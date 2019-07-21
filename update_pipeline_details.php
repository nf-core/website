<?php
//
// update_pipeline_details.php
// ---------------------------
// Pull the details of all nf-core pipelines from GitHub
// and write as JSON to a file.
// This script is called by deploy.php when the GitHub web hooks
// trigger due to an update.
//
// Note that the resulting file (public_html/pipelines.json) is
// ignored in the .gitignore file and will not be tracked in git history.
//
// Manual usage: on command line, simply execute this script:
//   $ php update_pipeline_details.php

// Load the twitter PHP library
require "includes/libraries/twitteroauth/autoload.php";
use Abraham\TwitterOAuth\TwitterOAuth;

// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);

// Get the twitter auth secrets
$config = parse_ini_file("config.ini");

// HTTP header to use on API GET requests
$api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => [
            'User-Agent: PHP',
            'Accept:application/vnd.github.mercy-preview+json' // Needed to get topics (keywords) for now
        ]
    ]
]);

// Final filename to write JSON to
$results_fn = dirname(__FILE__).'/public_html/pipelines.json';

// Load a copy of the existing JSON file, if it exists
$old_json = false;
$tweets = array();
if(file_exists($results_fn)){
    $old_json = json_decode(file_get_contents($results_fn), true);
}


// Function to sort assoc array by key value (name)
function sort_name($a,$b) {
    return strcmp($a["full_name"], $b["full_name"]);
}
// Function to sort assoc array by key value (datestamp)
function sort_datestamp($a,$b) {
    return strtotime($a['published_at']) - strtotime($b['published_at']);
}

// Initialise the results array with the current time and placeholders
$results = array(
    'updated' => time(),
    'pipeline_count' => 0,
    'published_count' => 0,
    'devel_count' => 0,
    'archived_count' => 0
);

// Fetch all repositories at nf-core
$gh_api_url = 'https://api.github.com/orgs/nf-core/repos?per_page=100';
$gh_repos = json_decode(file_get_contents($gh_api_url, false, $api_opts));
if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
    var_dump($http_response_header);
    die("Could not fetch nf-core repositories! $gh_api_url");
}

// Save data from non-ignored repositories
$ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
$ignored_topics = parse_ini_file("ignored_repos.ini")['topics'];
foreach($gh_repos as $repo){
    if(!in_array($repo->name, $ignored_repos)){
        $topics = array();
        if(!is_null($repo->topics)){
            foreach($repo->topics as $topic){
                if(!in_array($topic, $ignored_topics)){
                    $topics[] = $topic;
                }
            }
        }
        $results['remote_workflows'][] = array(
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
            'watchers_count' => $repo->watchers_count,
            'forks_count' => $repo->forks_count,
            'archived' => $repo->archived,
            'watchers' => $repo->watchers
        );
    }
}
// Sort workflows by name
usort($results['remote_workflows'], 'sort_name');

// Get additional release data for each repo
foreach($results['remote_workflows'] as $idx => $repo){
    // Fetch details of latest commits for this repo
    $results['remote_workflows'][$idx]['last_commit_sha'] = [ 'master' => '', 'dev' => '' ];
    $gh_master_refs_url = "https://api.github.com/repos/{$repo['full_name']}/git/refs/heads/master";
    $gh_master_commits = json_decode(file_get_contents($gh_master_refs_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        echo("Could not fetch nf-core master sha info! $gh_master_refs_url");
    }
    if(isset($gh_master_commits->object->sha)){
        $results['remote_workflows'][$idx]['last_commit_sha']['master'] = $gh_master_commits->object->sha;
    }
    $gh_dev_refs_url = "https://api.github.com/repos/{$repo['full_name']}/git/refs/heads/dev";
    $gh_dev_commits = json_decode(file_get_contents($gh_dev_refs_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        echo("Could not fetch nf-core dev sha info! $gh_dev_refs_url");
    }
    if(isset($gh_dev_commits->object->sha)){
        $results['remote_workflows'][$idx]['last_commit_sha']['dev'] = $gh_dev_commits->object->sha;
    }

    // Fetch release information for this repo
    $gh_releases_url = "https://api.github.com/repos/{$repo['full_name']}/releases";
    $gh_releases = json_decode(file_get_contents($gh_releases_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        echo("Could not fetch nf-core release info! $gh_releases_url");
    }

    // Save releases to results
    $results['remote_workflows'][$idx]['releases'] = [];
    foreach($gh_releases as $rel){
        $results['remote_workflows'][$idx]['releases'][] = array(
            'name' => $rel->name,
            'published_at' => $rel->published_at,
            'html_url' => $rel->html_url,
            'tag_name' => $rel->tag_name,
            'tag_sha' => NULL,
            'draft' => $rel->draft,
            'prerelease' => $rel->prerelease
        );
        if(strtotime($rel->published_at) > strtotime($results['remote_workflows'][$idx]['last_release'])){
            $results['remote_workflows'][$idx]['last_release'] = $rel->published_at;
        }
    }
    if(count($results['remote_workflows'][$idx]['releases']) > 0){
        // Sort releases by date, descending
        usort($results['remote_workflows'][$idx]['releases'], 'sort_datestamp');

        // Get commit hash information for each release
        $gh_tags_url = "https://api.github.com/repos/{$repo['full_name']}/tags";
        $gh_tags = json_decode(file_get_contents($gh_tags_url, false, $api_opts));
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            echo("Could not fetch nf-core tags info! $gh_tags_url");
        }
        foreach($gh_tags as $tag){
            foreach($results['remote_workflows'][$idx]['releases'] as $relidx => $rel){
                if($tag->name == $rel['tag_name']){
                    $results['remote_workflows'][$idx]['releases'][$relidx]['tag_sha'] = $tag->commit->sha;
                }
            }
        }
    }
}

// Count workflows
foreach($results['remote_workflows'] as $repo){
    $results['pipeline_count']++;
    if($repo['archived']){
        $results['archived_count']++;
    } else if(count($repo['releases']) > 0){
        $results['published_count']++;
    } else {
        $results['devel_count']++;
    }
}

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT)."\n";
file_put_contents($results_fn, $results_json);

////// Tweet about new releases
// Get old releases
$old_rel_tags = array();
foreach($old_json['remote_workflows'] as $old_pipeline){
    $old_rel_tags[$old_pipeline['name']] = array();
    // Collect releases from this pipeline
    foreach($old_pipeline['releases'] as $rel){
        if($rel['draft'] || $rel['prerelease']){
            continue;
        }
        $old_rel_tags[$old_pipeline['name']][] = $rel['tag_name'];
    }
}
// Go through new releases
foreach($results['remote_workflows'] as $new_pipeline){
    $rel_urls = array();
    foreach($new_pipeline['releases'] as $rel){
        if($rel['draft'] || $rel['prerelease']){
            continue;
        }
        // See if this tag name was in the previous JSON
        if(!in_array($rel['tag_name'], $old_rel_tags[$new_pipeline['name']])){
            // Prepare the tweet content!
            $tweet = 'Pipeline release! ';
            $tweet .= $new_pipeline['full_name'].' v'.$rel['tag_name'].' ('.$new_pipeline['description'].')';
            $tweet_url = "\n\nSee the changelog: ".$rel['html_url'];
            // 42 chars for tweet URL string with t.co url shortener
            while( (strlen($tweet) + 42) > $config['twitter_tweet_length'] ){
                $tweet = substr(rtrim($tweet, '.)'), 0, -3).'..)';
            }
            $tweets[] = $tweet.$tweet_url;
        }
    }
}

// Only tweet if we're on the live server!
if(count($tweets) > 0 && $_SERVER['SERVER_NAME'] == 'nf-co.re'){

    // Connect to twitter
    $connection = new TwitterOAuth(
        $config['twitter_key'],
        $config['twitter_secret'],
        $config['twitter_access_token'],
        $config['twitter_access_token_secret']
    );

    // Post the tweets
    foreach($tweets as $tweet){
        $post_tweets = $connection->post("statuses/update", ["status" => $tweet]);
    }
}

?>
