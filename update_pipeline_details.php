<?php
//
// update_pipeline_details.php
// ---------------------------
// Pull the details of all nf-core pipelines from GitHub
// and write as JSON to a file.
// This script is called by deploy_pipelines.php when the GitHub web hooks
// trigger due to an update.
//
// Note that the resulting file (public_html/pipelines.json) is
// ignored in the .gitignore file and will not be tracked in git history.
//
// Manual usage: on command line, simply execute this script:
//   $ php update_pipeline_details.php

echo "\n\nUpdating pipeline details - " . date('Y-m-d h:i:s') . "\n";

// Load the twitter PHP library
require 'vendor/autoload.php';

use Abraham\TwitterOAuth\TwitterOAuth;

// Allow PHP fopen to work with remote links
ini_set('allow_url_fopen', 1);

// Get the twitter auth secrets
$config = parse_ini_file('config.ini');

function get_gh_api($gh_api_url) {
    global $config;
    $gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

    // HTTP header to use on API GET requests
    $gh_api_opts = stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => [
                'User-Agent: PHP',
                'Accept:application/vnd.github.mercy-preview+json', // Needed to get topics (keywords) for now
                "Authorization: Basic $gh_auth",
            ],
        ],
    ]);

    // Get API response
    $first_page = true;
    $next_page = false;
    $counter = 1;
    $res = [];
    while ($first_page || $next_page) {
        // reset loop vars
        $first_page = false;
        // Get GitHub API results
        if ($next_page) {
            $gh_api_url = $next_page;
        }
        $tmp_results = json_decode(file_get_contents($gh_api_url, false, $gh_api_opts));
        if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
            die("\n-------- START ERROR " . date('Y-m-d h:i:s') . " --------\nCould not fetch $gh_api_url \n");
            var_dump($http_response_header);
            echo "\n$tmp_results\n";
            die("\n-------- END ERROR " . date('Y-m-d h:i:s') . " --------\nCould not fetch $gh_api_url \n");
        }

        array_push($res, ...$tmp_results);
        // Look for URL to next page of API results
        $m_array = preg_grep('/rel="next"/', $http_response_header);

        if (count($m_array) > 0) {
            $counter++;
            // check if page parameter is in gh_api_url
            if (strpos($gh_api_url, 'page=') === false) {
                $next_page = $gh_api_url . '?page=2';
            } else {
                $next_page = preg_replace('/page=\d+/', 'page=' . $counter, $gh_api_url);
            }
        } else {
            $next_page = false;
        }
    }
    return $res;
}

// Final filenames to write JSON to
$results_fn = dirname(__FILE__) . '/public_html/pipelines.json';
$pipeline_names_fn = dirname(__FILE__) . '/public_html/pipeline_names.json';

// Load a copy of the existing JSON file, if it exists
$old_json = false;
$tweets = [];
if (file_exists($results_fn)) {
    $old_json = json_decode(file_get_contents($results_fn), true);
}

// Function to sort assoc array by key value (name)
function sort_name($a, $b) {
    return strcmp($a['full_name'], $b['full_name']);
}
// Function to sort assoc array by key value (datestamp)
function sort_datestamp($a, $b) {
    return strtotime($a['published_at']) - strtotime($b['published_at']);
}

// Initialise the results array with the current time and placeholders
$results = [
    'updated' => time(),
    'pipeline_count' => 0,
    'published_count' => 0,
    'devel_count' => 0,
    'archived_count' => 0,
];

// Fetch all repositories at nf-core
$gh_repos = get_gh_api('https://api.github.com/orgs/nf-core/repos');

// Save data from non-ignored repositories
$ignored_repos = parse_ini_file('ignored_repos.ini')['repos'];
$ignored_topics = parse_ini_file('ignored_repos.ini')['topics'];
foreach ($gh_repos as $repo) {
    if (!in_array($repo->name, $ignored_repos)) {
        $topics = [];
        if (!is_null($repo->topics)) {
            foreach ($repo->topics as $topic) {
                if (!in_array($topic, $ignored_topics)) {
                    $topics[] = $topic;
                }
            }
        }
        $results['remote_workflows'][] = [
            'id' => $repo->id,
            'name' => $repo->name,
            'full_name' => $repo->full_name,
            'private' => $repo->private,
            'html_url' => $repo->html_url,
            'description' => $repo->description,
            'topics' => $topics,
            'created_at' => $repo->created_at,
            'updated_at' => $repo->updated_at,
            'pushed_at' => $repo->pushed_at,
            'last_release' => 0,
            'git_url' => $repo->git_url,
            'ssh_url' => $repo->ssh_url,
            'clone_url' => $repo->clone_url,
            'size' => $repo->size,
            'stargazers_count' => $repo->stargazers_count,
            'forks_count' => $repo->forks_count,
            'archived' => $repo->archived,
        ];
    }
}
// Sort workflows by name
usort($results['remote_workflows'], 'sort_name');

// Get additional release data for each repo
foreach ($results['remote_workflows'] as $idx => $repo) {
    // Fetch release information for this repo
    $gh_releases = get_gh_api("https://api.github.com/repos/{$repo['full_name']}/releases");

    // Save releases to results
    $results['remote_workflows'][$idx]['releases'] = [];
    foreach ($gh_releases as $rel) {
        // Skip if a draft release or prerelease
        if ($rel->draft) {
            continue;
        }
        if ($rel->prerelease) {
            continue;
        }

        $results['remote_workflows'][$idx]['releases'][] = [
            'name' => $rel->name,
            'published_at' => $rel->published_at,
            'html_url' => $rel->html_url,
            'tag_name' => $rel->tag_name,
            'tag_sha' => null,
            'draft' => $rel->draft,
            'prerelease' => $rel->prerelease,
            'tarball_url' => $rel->tarball_url,
            'zipball_url' => $rel->zipball_url,
        ];
        if (strtotime($rel->published_at) > strtotime($results['remote_workflows'][$idx]['last_release'])) {
            $results['remote_workflows'][$idx]['last_release'] = $rel->published_at;
        }
    }
    if (count($results['remote_workflows'][$idx]['releases']) > 0) {
        // Sort releases by date, descending
        usort($results['remote_workflows'][$idx]['releases'], 'sort_datestamp');

        // Get commit hash information for each release
        $gh_tags = get_gh_api("https://api.github.com/repos/{$repo['full_name']}/tags");
        foreach ($gh_tags as $tag) {
            foreach ($results['remote_workflows'][$idx]['releases'] as $relidx => $rel) {
                if ($tag->name == $rel['tag_name']) {
                    $results['remote_workflows'][$idx]['releases'][$relidx]['tag_sha'] = $tag->commit->sha;
                }
            }
        }
    }
}

// Count workflows
$pipeline_names = [];
foreach ($results['remote_workflows'] as $repo) {
    $results['pipeline_count']++;
    if ($repo['archived']) {
        $results['archived_count']++;
    } elseif (count($repo['releases']) > 0) {
        $results['published_count']++;
        $pipeline_names[] = $repo['name'];
    } else {
        $results['devel_count']++;
        $pipeline_names[] = $repo['name'];
    }
}

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT) . "\n";
file_put_contents($results_fn, $results_json);

// Print simple list of pipelines to a file
file_put_contents($pipeline_names_fn, json_encode(['pipeline' => $pipeline_names]));

////// Tweet about new releases
// Get old releases
$old_rel_tags = [];
if ($old_json) {
    foreach ($old_json['remote_workflows'] as $old_pipeline) {
        $old_rel_tags[$old_pipeline['name']] = [];
        // Collect releases from this pipeline
        foreach ($old_pipeline['releases'] as $rel) {
            if ($rel['draft'] || $rel['prerelease']) {
                continue;
            }
            $old_rel_tags[$old_pipeline['name']][] = $rel['tag_name'];
        }
    }
}
// Go through new releases
foreach ($results['remote_workflows'] as $new_pipeline) {
    $rel_urls = [];
    foreach ($new_pipeline['releases'] as $rel) {
        if ($rel['draft'] || $rel['prerelease']) {
            continue;
        }
        // See if this tag name was in the previous JSON
        if (
            !isset($old_rel_tags[$new_pipeline['name']]) ||
            !in_array($rel['tag_name'], $old_rel_tags[$new_pipeline['name']])
        ) {
            // Prepare the tweet content!
            $tweet = 'Pipeline release! ';
            $tweet .= $new_pipeline['full_name'] . ' v' . $rel['tag_name'] . ' (' . $new_pipeline['description'] . ')';
            $tweet_url = "\n\nSee the changelog: " . $rel['html_url'];
            // 42 chars for tweet URL string with t.co url shortener
            while (strlen($tweet) + 42 > $config['twitter_tweet_length']) {
                $tweet = substr(rtrim($tweet, '.)'), 0, -3) . '..)';
            }
            $tweets[] = $tweet . $tweet_url;
            echo 'Found new release for ' . $new_pipeline['full_name'] . ': ' . $rel['tag_name'] . "\n";
        }
    }
}

// Only tweet if we're on the live server!
if (count($tweets) > 0) {
    if (isset($config['twitter_key']) && $config['twitter_key'] != 'TWITTER_KEY') {
        // Connect to twitter
        $connection = new TwitterOAuth(
            $config['twitter_key'],
            $config['twitter_secret'],
            $config['twitter_access_token'],
            $config['twitter_access_token_secret'],
        );
        // loop through pages of twitter api results
        $page = 0;
        $max_pages = 100;

        $already_tweeted = false;
        foreach ($tweets as $tweet) {
            while ($page < $max_pages) {
                $page++;
                $params = [
                    'count' => 200,
                    'page' => $page,
                    'exclude_replies' => true,
                    'exclude_retweets' => true,
                ];
                $old_tweets = $connection->get('statuses/user_timeline', $params);
                // Check for substring of the new tweet in old tweets
                foreach ($old_tweets as $old_tweet) {
                    preg_match('/Pipeline release!.*?\(/', $tweet, $match);
                    if (strpos($old_tweet->text, $match[0]) !== false) {
                        $already_tweeted = true;
                        echo 'Already tweeted at https://twitter.com/nf_core/status/' .
                            $old_tweet->id .
                            "! Not sending it.\n";
                        break;
                    }
                }

                if (count($old_tweets) == 0 || $already_tweeted) {
                    break;
                }
            }
            if (!$already_tweeted) {
                $connection->post('statuses/update', ['status' => $tweet]);
                echo "Sent tweet: $tweet\n";
            }
        }
    } else {
        echo "Not sending tweets because config twitter_key is not set.\n";
    }
}

echo "\nupdate_pipeline_details done " . date('Y-m-d h:i:s') . "\n\n";
