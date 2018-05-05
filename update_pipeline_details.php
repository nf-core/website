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


// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);
// HTTP header to use on API GET requests
$api_opts = stream_context_create(['http' => ['method' => 'GET', 'header' => ['User-Agent: PHP']]]);
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
    var_dump($gh_repos);
    die("Could not fetch nf-core repositories! $gh_api_url");
}

// Save data from non-ignored repositories
$ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
foreach($gh_repos as $repo){
    if(!in_array($repo->name, $ignored_repos)){
        $results['remote_workflows'][] = array(
            'id' => $repo->id,
            'name' => $repo->name,
            'full_name' => $repo->full_name,
            'private' => $repo->private,
            'html_url' => $repo->html_url,
            'description' => $repo->description,
            'created_at' => $repo->created_at,
            'updated_at' => $repo->updated_at,
            'pushed_at' => $repo->pushed_at,
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
    // Fetch release information for this repo
    $gh_releases_url = "https://api.github.com/repos/{$repo['full_name']}/releases";
    $gh_releases = json_decode(file_get_contents($gh_releases_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($gh_releases);
        die("Could not fetch nf-core release info! $gh_releases_url");
    }

    // Save releases to results
    $results['remote_workflows'][$idx]['releases'] = [];
    foreach($gh_releases as $rel){
        $results['remote_workflows'][$idx]['releases'][] = array(
            'name' => $rel->name,
            'published_at' => $rel->published_at,
            'html_url' => $rel->html_url,
            'tag_name' => $rel->tag_name,
            'tag_sha' => None,
            'draft' => $rel->draft,
            'prerelease' => $rel->prerelease
        );
    }
    if(count($results['remote_workflows'][$idx]['releases']) > 0){
        // Sort releases by date, descending
        usort($results['remote_workflows'][$idx]['releases'], 'sort_datestamp');

        // Get commit hash information for each release
        $gh_tags_url = "https://api.github.com/repos/{$repo['full_name']}/tags";
        $gh_tags = json_decode(file_get_contents($gh_tags_url, false, $api_opts));
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($gh_tags);
            die("Could not fetch nf-core tags info! $gh_tags_url");
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
$results_fn = dirname(__FILE__).'/public_html/pipelines.json';
file_put_contents($results_fn, $results_json)

?>
