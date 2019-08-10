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

$debug = false;
$num_api_calls = 0;
$max_repos = false;
$max_comments = false;

// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);

// Use same updated time for everything
$updated = time();


// Final filename to write JSON to
$results_fn = dirname(__FILE__).'/nfcore_issue_stats.json';
$contribs_fn_root = dirname(__FILE__).'/contributor_stats/';

// Load a copy of the existing JSON file, if it exists
if(file_exists($results_fn)){
    if($debug){ echo "Loading previous stats: $results_fn\n"; }
    $results = json_decode(file_get_contents($results_fn), true);
}
$results['updated'] = $updated;
$base_stats = [
    'count' => 0,
    'open_count' => 0,
    'comments_count' => 0,
    'authors_count' => [],
    'median_close_time' => 0,
    'median_response_time' => 0
];
$results['stats'][$updated]['issues'] = $base_stats;
$results['stats'][$updated]['prs'] = $base_stats;
$results['stats']['issues']['daily_opened'] = [];
$results['stats']['issues']['daily_closed'] = [];
$results['stats']['prs']['daily_opened'] = [];
$results['stats']['prs']['daily_closed'] = [];

// Get auth secrets
$config = parse_ini_file("config.ini");
$gh_auth = base64_encode($config['github_username'].':'.$config['github_access_token']);


//
//
// GitHub API calls
//
//


// HTTP header to use on GitHub API GET requests
$gh_api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => [
            'User-Agent: PHP',
            "Authorization: Basic $gh_auth"
        ]
    ]
]);

// Load details of the pipelines
$pipelines_json = json_decode(file_get_contents('public_html/pipelines.json'), true);
$pipelines = $pipelines_json['remote_workflows'];
// Get list of other repos
$repos = parse_ini_file("ignored_repos.ini")['repos'];
foreach($pipelines as $pipeline){
    $repos[] = $pipeline['name'];
}


// Get all issues for all repos

// Get the current number of organisation members
// Remember pagination!
$num_repos = 0;
foreach($repos as $repo){
    $num_repos += 1;
    if(is_numeric($max_repos)){
        if($num_repos > $max_repos){
            break;
        }
    }
    if(!isset($results['repos'][$repo])){
        $results['repos'][$repo] = [
            'issues' => [],
            'prs' => [],
        ];
    }
    $gh_issues_url = 'https://api.github.com/repos/nf-core/'.$repo.'/issues?state=all';
    $first_page = true;
    $next_page = false;
    while($first_page || $next_page){
        // reset loop vars
        $first_page = false;
        // Get GitHub API results
        if($next_page){
            $gh_issues_url = $next_page;
        }
        if($debug){ echo "Fetching $gh_issues_url\n"; }
        $num_api_calls += 1;
        $gh_issues = json_decode(file_get_contents($gh_issues_url, false, $gh_api_opts), true);
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core/$repo issues! $gh_issues_url");
        }
        // Don't clobber what we already have - set keys individually
        foreach($gh_issues as $issue){
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
            if($issue['closed_at']){
                $closed_wait = strtotime($issue['closed_at']) - strtotime($issue['created_at']);
                $results['repos'][$repo][$issue_type][$id]['closed_wait'] = $closed_wait;
                $results['stats'][$issue_type]['close_times'][] = $closed_wait;
            }
            $results['repos'][$repo][$issue_type][$id]['created_by'] = $author;
            // Contributor stats
            if(!isset($results['authors'][$author][$issue_type]['num_created'])){
                $results['authors'][$author][$issue_type]['num_created'] = 0;
            }
            $results['authors'][$author][$issue_type]['num_created'] += 1;
            // Global stats
            $results['stats'][$updated][$issue_type]['count'] += 1;
            if($issue['state'] == 'open'){
                $results['stats'][$updated][$issue_type]['open_count'] += 1;
            }
            $results['stats'][$updated][$issue_type]['comments_count'] += $issue['comments'];
            $results['stats'][$updated][$issue_type]['authors_count'][$author] = 0;
            // Daily opened / closed stats
            $created_at_day = date('d-m-Y', strtotime($issue['created_at']));
            if(!isset($results['stats'][$issue_type]['daily_opened'][$created_at_day])){
                $results['stats'][$issue_type]['daily_opened'][$created_at_day] = 0;
            }
            $results['stats'][$issue_type]['daily_opened'][$created_at_day] += 1;
            if($issue['closed_at']){
                $closed_at_day = date('d-m-Y', strtotime($issue['closed_at']));
                if(!isset($results['stats'][$issue_type]['daily_closed'][$closed_at_day])){
                    $results['stats'][$issue_type]['daily_closed'][$closed_at_day] = 0;
                }
                $results['stats'][$issue_type]['daily_closed'][$created_at_day] += 1;
            }
        }
        // Look for URL to next page of API results
        $next_page = false;
        $m_array = preg_grep('/rel="next"/', $http_response_header);
        if(count($m_array) > 0){
            preg_match('/<([^>]+)>; rel="next"/', array_values($m_array)[0], $matches);
            if(isset($matches[1])){
                $next_page = $matches[1];
            }
        }
    }
}



// Get details of first comment for every issue
$num_calls = 0;
foreach($repos as $repo){
    foreach(['prs', 'issues'] as $issue_type){
        if(!isset($results['repos'][$repo][$issue_type])){
            continue;
        }
        foreach($results['repos'][$repo][$issue_type] as $id => $issue){
            // For developing - if $max_calls is set, check if we should bail
            $num_calls += 1;
            if(is_numeric($max_comments)){
                if($num_calls > $max_comments){
                    break 2; // break the outer loop too
                }
            }
            // Skip issues with no comments
            if($issue['num_comments'] == 0){
                continue;
            }
            // Lots of API calls! Only look for this if we don't already have details
            if(isset($issue['first_reply'])){
                continue;
            }
            // Some issues only have comments from the author
            if(isset($issue['talking_to_myself']) && $issue['talking_to_myself'] == $issue['num_comments']){
                continue;
            }
            // Ok, should just be issues with new comments from here on
            $gh_comments_url = $issue['comments_url'];
            if($debug){ echo "Fetching $gh_comments_url\n"; }
            $num_api_calls += 1;
            $gh_comments = json_decode(file_get_contents($gh_comments_url, false, $gh_api_opts), true);
            if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
                var_dump($http_response_header);
                die("Could not fetch nf-core/$repo issue #$id! $gh_comments_url");
            }
            $comments_by_created = [];
            $comment_timestamps = [];
            foreach($gh_comments as $comment){
                // Only save first comment if it wasn't the issue author
                $ts = strtotime($comment['created_at']);
                $author = $comment['user']['login'];
                if($issue['created_by'] !== $author){
                    $comment_timestamps[] = $ts;
                    $comments_by_created[$ts] = $comment;
                }
                // Contributor stats
                if(!isset($results['authors'][$author][$issue_type]['num_replies'])){
                    $results['authors'][$author][$issue_type]['num_replies'] = 0;
                }
                $results['authors'][$author][$issue_type]['num_replies'] += 1;
                // Global stats
                $results['stats'][$updated][$issue_type]['authors_count'][$author] = 0;
            }
            // Special case - the only comments are from the original issue author
            if(count($comment_timestamps) == 0){
                $results['repos'][$repo][$issue_type][$id]['talking_to_myself'] = $issue['num_comments'];
            } else {
                $first_comment_ts = min($comment_timestamps);
                $results['repos'][$repo][$issue_type][$id]['first_reply'] = $comments_by_created[$first_comment_ts]['created_at'];
                $response_wait = $first_comment_ts - strtotime($issue['created_at']);
                $results['repos'][$repo][$issue_type][$id]['first_reply_wait'] = $response_wait;
                $results['stats'][$issue_type]['response_times'][] = $response_wait;
                // Contributor - first responder
                $first_responder = $comments_by_created[$first_comment_ts]['user']['login'];
                $results['repos'][$repo][$issue_type][$id]['first_reply_by'] = $first_responder;
                if(!isset($results['authors'][$first_responder][$issue_type]['num_first_response'])){
                    $results['authors'][$first_responder][$issue_type]['num_first_response'] = 0;
                }
                $results['authors'][$first_responder][$issue_type]['num_first_response'] += 1;
            }
        }
    }
}

// Count the author stats
$results['stats'][$updated]['issues']['authors_count'] = count($results['stats'][$updated]['issues']['authors_count']);
$results['stats'][$updated]['prs']['authors_count'] = count($results['stats'][$updated]['prs']['authors_count']);

// Calculate the median times
function array_median($arr){
    arsort($arr);
    $keys = array_keys($arr);
    return $arr[$keys[floor(count($keys)/2)]];
}
$results['stats'][$updated]['issues']['median_close_time'] = array_median($results['stats']['issues']['close_times']);
$results['stats'][$updated]['issues']['median_response_time'] = array_median($results['stats']['issues']['response_times']);
$results['stats'][$updated]['prs']['median_close_time'] = array_median($results['stats']['prs']['close_times']);
$results['stats'][$updated]['prs']['median_response_time'] = array_median($results['stats']['prs']['response_times']);

//
//
// DONE - save results to JSON file
//
//

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT)."\n";
file_put_contents($results_fn, $results_json);

echo("\nupdate_issue_stats done with $num_api_calls API calls - ".date("Y-m-d h:i:s")."\n\n");
