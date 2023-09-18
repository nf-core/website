<?php
//
// nfcore_issue_stats.json
// ---------------------------
// To get nice statistics like how fast we respond to issues and who
// creates and responds to issues, we query the GitHub API to get issue stats.
//
// Note that the resulting file (nfcore_issue_stats.json) is
// ignored in the .gitignore file and will not be tracked in git history.
//
// Manual usage: on command line, simply execute this script:
//   $ php update_issue_stats.php

$debug = true;
$num_api_calls = 0;
$max_repos = false;
$max_comments = false;

// Allow PHP fopen to work with remote links
ini_set('allow_url_fopen', 1);

// Use same updated time for everything
$updated = time();

// Final filename to write JSON to
$results_fn = dirname(__FILE__) . '/nfcore_issue_stats.json';
$results = ['updated' => $updated];
// Load a copy of the existing JSON file, if it exists
if (file_exists($results_fn)) {
    if ($debug) {
        echo "Loading previous stats: $results_fn\n";
    }
    $results = json_decode(file_get_contents($results_fn), true);
}
$prev_updated = $results['updated'];
$results['updated'] = $updated;

// First check if we've pulled these in the past hour
$skip_issue_update = false;
// if ($updated - $prev_updated < 60 * 60) {
//     $results['updated'] = $prev_updated;
//     $updated = $prev_updated;
//     $skip_issue_update = true;
//     if ($debug) {
//         echo "Skipping repo comments pull as cache is from past hour\n";
//     }
// } elseif ($debug) {
//     echo 'Issues results are ' . ($updated - $prev_updated) . " seconds old - pulling again\n";
// }

if (!$skip_issue_update) {
    $base_stats = [
        'count' => 0,
        'open_count' => 0,
        'comments_count' => 0,
        'authors_count' => [],
        'median_close_time' => 0,
        'median_response_time' => 0,
    ];
    $results['stats'][$updated]['issues'] = $base_stats;
    $results['stats'][$updated]['prs'] = $base_stats;
    $results['stats']['issues']['daily_opened'] = [];
    $results['stats']['issues']['daily_closed'] = [];
    $results['stats']['issues']['close_times'] = [];
    $results['stats']['prs']['daily_opened'] = [];
    $results['stats']['prs']['daily_closed'] = [];
    $results['stats']['prs']['close_times'] = [];
}

// Get auth secrets
$config = parse_ini_file('config.ini');
$gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

//
//
// GitHub API calls
//
//

// HTTP header to use on GitHub API GET requests
$gh_api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => ['User-Agent: PHP', "Authorization: Basic $gh_auth"],
    ],
]);

// Load details of the pipelines
$pipelines_json = json_decode(file_get_contents(dirname(__FILE__) . '/public_html/pipelines.json'), true);
$pipelines = $pipelines_json['remote_workflows'];
// Get list of other repos
$repos = []; //parse_ini_file('ignored_repos.ini')['repos'];
foreach ($pipelines as $pipeline) {
    $repos[] = $pipeline['name'];
}

// Delete cached pipelines stats for pipelines that have been deleted
if (array_key_exists('repos', $results)) {
    foreach (array_keys($results['repos']) as $repo_name) {
        if (!in_array($repo_name, $repos)) {
            echo "\nRemoving $repo_name from the cached results as it appears to have been deleted.\n";
            unset($results['repos'][$repo_name]);
        }
    }
}

// Get all issues for all repos

// Get the current number of organisation members
// Remember pagination!
$num_repos = 0;
foreach ($repos as $repo) {
    $num_repos += 1;
    // First check if we've pulled these in the past hour and skip if so
    if ($skip_issue_update) {
        break;
    }

    if (is_numeric($max_repos)) {
        if ($num_repos > $max_repos) {
            break;
        }
    }
    if (!isset($results['repos'][$repo])) {
        $results['repos'][$repo] = [
            'issues' => [],
            'prs' => [],
        ];
    }
    $gh_issues_url = 'https://api.github.com/repos/sanger-tol/' . $repo . '/issues?state=all';
    $first_page = true;
    $next_page = false;
    while ($first_page || $next_page) {
        // reset loop vars
        $first_page = false;
        // Get GitHub API results
        if ($next_page) {
            $gh_issues_url = $next_page;
        }
        if ($debug) {
            echo "Fetching $gh_issues_url\n";
        }
        $num_api_calls += 1;
        $gh_issues = json_decode(file_get_contents($gh_issues_url, false, $gh_api_opts), true);
        if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
            var_dump($http_response_header);
            echo "\nCould not fetch nf-core/$repo issues! $gh_issues_url\n";
            continue;
        }
        // Don't clobber what we already have - set keys individually
        foreach ($gh_issues as $issue) {
            // Set up array keys
            $issue_type = isset($issue['pull_request']) ? 'prs' : 'issues';
            $id = $issue['number'];
            $author = $issue['user']['login'];
            // Save fields of interest
            $results['repos'][$repo][$issue_type][$id]['url'] = $issue['url'];
            $results['repos'][$repo][$issue_type][$id]['comments_url'] = $issue['comments_url'];
            $results['repos'][$repo][$issue_type][$id]['html_url'] = $issue['html_url'];
            $results['repos'][$repo][$issue_type][$id]['state'] = $issue['state'];
            $results['repos'][$repo][$issue_type][$id]['num_comments'] = $issue['comments'];
            $results['repos'][$repo][$issue_type][$id]['created_at'] = $issue['created_at'];
            $results['repos'][$repo][$issue_type][$id]['updated_at'] = $issue['updated_at'];
            $results['repos'][$repo][$issue_type][$id]['closed_at'] = $issue['closed_at'];
            if ($issue['closed_at']) {
                $closed_wait = strtotime($issue['closed_at']) - strtotime($issue['created_at']);
                $results['repos'][$repo][$issue_type][$id]['closed_wait'] = $closed_wait;
                $results['stats'][$issue_type]['close_times'][] = $closed_wait;
            }
            $results['repos'][$repo][$issue_type][$id]['created_by'] = $author;
            // Contributor stats
            if (!isset($results['authors'][$author][$issue_type]['num_created'])) {
                $results['authors'][$author][$issue_type]['num_created'] = 0;
            }
            $results['authors'][$author][$issue_type]['num_created'] += 1;
            if (!isset($results['authors'][$author]['first_contribution'])) {
                $results['authors'][$author]['first_contribution'] = strtotime($issue['created_at']);
            }
            $results['authors'][$author]['first_contribution'] = min(
                $results['authors'][$author]['first_contribution'],
                strtotime($issue['created_at']),
            );
            // Global stats
            $results['stats'][$updated][$issue_type]['count'] += 1;
            if ($issue['state'] == 'open') {
                $results['stats'][$updated][$issue_type]['open_count'] += 1;
            }
            $results['stats'][$updated][$issue_type]['comments_count'] += $issue['comments'];
            $results['stats'][$updated][$issue_type]['authors'][$author] = 0;
            // Daily opened / closed stats
            $created_at_day = date('d-m-Y', strtotime($issue['created_at']));
            if (!isset($results['stats'][$issue_type]['daily_opened'][$created_at_day])) {
                $results['stats'][$issue_type]['daily_opened'][$created_at_day] = 0;
            }
            $results['stats'][$issue_type]['daily_opened'][$created_at_day] += 1;
            if ($issue['closed_at']) {
                $closed_at_day = date('d-m-Y', strtotime($issue['closed_at']));
                if (!isset($results['stats'][$issue_type]['daily_closed'][$closed_at_day])) {
                    $results['stats'][$issue_type]['daily_closed'][$closed_at_day] = 0;
                }
                $results['stats'][$issue_type]['daily_closed'][$created_at_day] += 1;
            }
        }
        // Look for URL to next page of API results
        $next_page = false;
        $m_array = preg_grep('/rel="next"/', $http_response_header);
        if (count($m_array) > 0) {
            preg_match('/<([^>]+)>; rel="next"/', array_values($m_array)[0], $matches);
            if (isset($matches[1])) {
                $next_page = $matches[1];
            }
        }
    }
}

// Get details of first comment for every issue
$num_calls = 0;
foreach ($repos as $repo) {
    foreach (['prs', 'issues'] as $issue_type) {
        if (!isset($results['repos'][$repo][$issue_type])) {
            continue;
        }
        foreach ($results['repos'][$repo][$issue_type] as $id => $issue) {
            // For developing - if $max_calls is set, check if we should bail
            $num_calls += 1;
            if (is_numeric($max_comments)) {
                if ($num_calls > $max_comments) {
                    break 2; // break the outer loop too
                }
            }

            // Skip issues with no comments
            if ($issue['num_comments'] == 0) {
                continue;
            }

            // Check we if we have pulled this already
            $comments_fn = dirname(__FILE__) . '/api_cache/issue_comments/' . $repo . '/' . $id . '.json';
            $fetch_json = true;
            if (file_exists($comments_fn)) {
                if ($debug) {
                    echo "Loading previous comments: $comments_fn\n";
                }
                $gh_comments = json_decode(file_get_contents($comments_fn), true);
                // Check that there haven't been new comments since we cached this
                if ($issue['num_comments'] == count($gh_comments)) {
                    $fetch_json = false;
                } elseif ($debug) {
                    echo "New comments for nf-core/$repo issue #$id: now " .
                        htmlspecialchars($issue['num_comments'], ENT_QUOTES, 'UTF-8') .
                        ', previously had ' .
                        count($gh_comments) .
                        "\n";
                }
            }
            if ($fetch_json) {
                $gh_comments = [];
                $gh_comments_url = $issue['comments_url'];
                $first_page = true;
                $next_page = false;
                while ($first_page || $next_page) {
                    $first_page = false;
                    if ($next_page) {
                        $gh_comments_url = $next_page;
                    }
                    if ($debug) {
                        echo "Fetching $gh_comments_url\n";
                    }
                    $num_api_calls += 1;
                    $gh_new_comments = json_decode(file_get_contents($gh_comments_url, false, $gh_api_opts), true);
                    $gh_comments = array_merge($gh_comments, $gh_new_comments);
                    if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
                        var_dump($http_response_header);
                        echo "\nCould not fetch nf-core/$repo issue #$id! $gh_comments_url\n";
                        continue;
                    }
                    // Look for URL to next page of API results
                    $next_page = false;
                    $m_array = preg_grep('/rel="next"/', $http_response_header);
                    if (count($m_array) > 0) {
                        preg_match('/<([^>]+)>; rel="next"/', array_values($m_array)[0], $matches);
                        if (isset($matches[1])) {
                            $next_page = $matches[1];
                        }
                    }
                }
                // Save for next time
                if (!file_exists(dirname($comments_fn))) {
                    mkdir(dirname($comments_fn), 0777, true);
                }
                $gh_comments_json = json_encode($gh_comments, JSON_PRETTY_PRINT) . "\n";
                file_put_contents($comments_fn, $gh_comments_json);
            }

            // Ok, lets go through the comments
            $comments_by_created = [];
            $comment_timestamps = [];
            foreach ($gh_comments as $comment) {
                // Only save first comment if it wasn't the issue author
                $ts = strtotime($comment['created_at']);
                $author = $comment['user']['login'];
                if ($issue['created_by'] !== $author) {
                    $comment_timestamps[] = $ts;
                    $comments_by_created[$ts] = $comment;
                }
                // Contributor stats
                if (!isset($results['authors'][$author][$issue_type]['num_replies'])) {
                    $results['authors'][$author][$issue_type]['num_replies'] = 0;
                }
                $results['authors'][$author][$issue_type]['num_replies'] += 1;
                if (!isset($results['authors'][$author]['first_contribution'])) {
                    $results['authors'][$author]['first_contribution'] = strtotime($issue['created_at']);
                }
                $results['authors'][$author]['first_contribution'] = min(
                    $results['authors'][$author]['first_contribution'],
                    strtotime($issue['created_at']),
                );
                // Global stats
                $results['stats'][$updated][$issue_type]['authors'][$author] = 0;
            }
            // Special case - the only comments are from the original issue author
            if (count($comment_timestamps) == 0) {
                $results['repos'][$repo][$issue_type][$id]['talking_to_myself'] = $issue['num_comments'];
            } else {
                $first_comment_ts = min($comment_timestamps);
                $results['repos'][$repo][$issue_type][$id]['first_reply'] =
                    $comments_by_created[$first_comment_ts]['created_at'];
                $response_wait = $first_comment_ts - strtotime($issue['created_at']);
                $results['repos'][$repo][$issue_type][$id]['first_reply_wait'] = $response_wait;
                $results['stats'][$issue_type]['response_times'][] = $response_wait;
                // Contributor - first responder
                $first_responder = $comments_by_created[$first_comment_ts]['user']['login'];
                $results['repos'][$repo][$issue_type][$id]['first_reply_by'] = $first_responder;
                if (!isset($results['authors'][$first_responder][$issue_type]['num_first_response'])) {
                    $results['authors'][$first_responder][$issue_type]['num_first_response'] = 0;
                }
                $results['authors'][$first_responder][$issue_type]['num_first_response'] += 1;
            }
        }
    }
}

// Count the author stats
// NB: This will be a bit wrong if we haven't grabbed issues, as won't be counting people who created issues but haven't had a comment
if (array_key_exists('stats', $results)) {
    $results['stats'][$updated]['issues']['authors_count'] = count($results['stats'][$updated]['issues']['authors']);
    $results['stats'][$updated]['prs']['authors_count'] = count($results['stats'][$updated]['prs']['authors']);
    unset($results['stats'][$updated]['issues']['authors']);
    unset($results['stats'][$updated]['prs']['authors']);

    // Calculate the median times
    function array_median($arr) {
        arsort($arr);
        $keys = array_keys($arr);
        return $arr[$keys[floor(count($keys) / 2)]];
    }
    $results['stats'][$updated]['issues']['median_close_time'] = array_median(
        $results['stats']['issues']['close_times'],
    );
    $results['stats'][$updated]['issues']['median_response_time'] = array_median(
        $results['stats']['issues']['response_times'],
    );
    $results['stats'][$updated]['prs']['median_close_time'] = array_median($results['stats']['prs']['close_times']);
    $results['stats'][$updated]['prs']['median_response_time'] = array_median(
        $results['stats']['prs']['response_times'],
    );
}
//
//
// DONE - save results to JSON file
//
//

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT) . "\n";
file_put_contents($results_fn, $results_json);

echo "\nupdate_issue_stats done with $num_api_calls API calls - " . date('Y-m-d h:i:s') . "\n\n";
