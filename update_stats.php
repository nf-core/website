<?php
//
// nfcore_stats.json
// ---------------------------
// GitHub shows traffic in the form of repo views and clones, however
// the data is only available for two weeks.
// We want it forever! So this script scrapes and saves the data.
// It is intended to be run routinely using a cronjob
//
// Note that the resulting file (nfcore_stats.json) is
// ignored in the .gitignore file and will not be tracked in git history.
//
// Manual usage: on command line, simply execute this script:
//   $ php nfcore_stats.json


// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);

// Get the GitHub auth secrets
$config = parse_ini_file("config.ini");
$auth = base64_encode($config['github_username'].':'.$config['github_access_token']);

// HTTP header to use on API GET requests
$api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => [
            'User-Agent: PHP',
            "Authorization: Basic $auth"
        ]
    ]
]);

// Final filename to write JSON to
$results_fn = dirname(__FILE__).'/nfcore_stats.json';

// Initialise the results array with the current time and placeholders
$results = array(
    'updated' => time(),
    'pipelines' => array()
);

// Load a copy of the existing JSON file, if it exists
if(file_exists($results_fn)){
    $results = json_decode(file_get_contents($results_fn), true);
}

// Load details of the pipelines
$pipelines_json = json_decode(file_get_contents('public_html/pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
$contribs_try_again = [];

// Fetch new statistics for each repo
foreach($pipelines as $wf){
    // Views
    $gh_views_url = 'https://api.github.com/repos/'.$wf->full_name.'/traffic/views';
    $gh_views = json_decode(file_get_contents($gh_views_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        die("Could not fetch nf-core repo views! $gh_views_url");
    }
    foreach($gh_views->views as $view){
        $results['pipelines'][$wf->name]['views_count'][$view->timestamp] = $view->count;
        $results['pipelines'][$wf->name]['views_uniques'][$view->timestamp] = $view->uniques;
    }
    // Clones
    $gh_clones_url = 'https://api.github.com/repos/'.$wf->full_name.'/traffic/clones';
    $gh_clones = json_decode(file_get_contents($gh_clones_url, false, $api_opts));
    if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        die("Could not fetch nf-core repo clones! $gh_clones_url");
    }
    foreach($gh_clones->clones as $clone){
        $results['pipelines'][$wf->name]['clones_count'][$clone->timestamp] = $clone->count;
        $results['pipelines'][$wf->name]['clones_uniques'][$clone->timestamp] = $clone->uniques;
    }
    // Contributors
    $gh_contributors_url = 'https://api.github.com/repos/'.$wf->full_name.'/stats/contributors';
    $gh_contributors = json_decode(file_get_contents($gh_contributors_url, false, $api_opts));
    // If the data hasn't been cached when you query a repository's statistics, you'll receive a 202 response;
    // a background job is also fired to start compiling these statistics.
    // Give the job a few moments to complete, and then submit the request again
    if(in_array("HTTP/1.1 202 Accepted", $http_response_header)){
        $contribs_try_again[$wf->name] = $gh_contributors_url;
    } else if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
        var_dump($http_response_header);
        die("Could not fetch nf-core repo contributors! $gh_contributors_url");
    }
    $results['pipelines'][$wf->name]['contributors'] = $gh_contributors;
    $results['pipelines'][$wf->name]['num_contributors'] = count($gh_contributors);

    // Recalculate totals
    foreach(['views_count', 'views_uniques', 'clones_count', 'clones_uniques'] as $ctype){
        $results['pipelines'][$wf->name][$ctype.'_total'] = 0;
        if(count($results['pipelines'][$wf->name][$ctype]) > 0){
            foreach($results['pipelines'][$wf->name][$ctype] as $stat){
                $results['pipelines'][$wf->name][$ctype.'_total'] += $stat;
            }
        }
    }
}

// Try contribs again now that we've let it fire
if(count($contribs_try_again) > 0){
    sleep(10);
    foreach($contribs_try_again as $wfname => $gh_contributors_url){
        $gh_contributors = json_decode(file_get_contents($gh_contributors_url, false, $api_opts));
        if(in_array("HTTP/1.1 202 Accepted", $http_response_header)){
            echo("Tried getting contributors after delay for $wfname, but took too long.");
        } else if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core repo contributors! $gh_contributors_url");
        }
        $results['pipelines'][$wfname]['contributors'] = $gh_contributors;
        $results['pipelines'][$wfname]['num_contributors'] = count($gh_contributors);
    }
}

// The data for commits per week is massive - remove it
foreach($pipelines as $wf){
    foreach($results['pipelines'][$wf->name]['contributors'] as $idx => $contributor){
        unset($results['pipelines'][$wf->name]['contributors'][$idx]->weeks);
    }
}

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT)."\n";
file_put_contents($results_fn, $results_json);

echo("update_pipeline_stats done " . mktime());
