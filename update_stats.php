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

// Use same updated time for everything
$updated = time();


// Final filename to write JSON to
$results_fn = dirname(__FILE__).'/nfcore_stats.json';

// Initialise the results array with the current time and placeholders
$results = array(
    'updated' => $updated,
    'pipelines' => array(),
    'core_repos' => array(),
    'slack' => array(),
    'gh_org_members' => array()
);

// Load a copy of the existing JSON file, if it exists
if(file_exists($results_fn)){
    $results = json_decode(file_get_contents($results_fn), true);
}
$results['updated'] = $updated;

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
$pipelines_json = json_decode(file_get_contents('public_html/pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
$contribs_try_again = [];

// Build array of repos to query
foreach($pipelines as $wf){
    if(!isset($results['pipelines'][$wf->name])){
        $results['pipelines'][$wf->name] = array();
    }
}
$ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
foreach($ignored_repos as $name){
    if(!isset($results['core_repos'][$name])){
        $results['core_repos'][$name] = array();
    }
}

// Get snapshot of key metrics for all repos

// Get the current number of organisation members
$gh_members_url = 'https://api.github.com/orgs/nf-core/members';
$gh_members = json_decode(file_get_contents($gh_members_url, false, $gh_api_opts));
if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
    var_dump($http_response_header);
    die("Could not fetch nf-core members! $gh_members_url");
}
$results['gh_org_members'][$updated] = count($gh_members);

// Fetch all repositories at nf-core
$gh_repos_url = 'https://api.github.com/orgs/nf-core/repos?per_page=100';
$gh_repos = json_decode(file_get_contents($gh_repos_url, false, $gh_api_opts));
if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
    var_dump($http_response_header);
    die("Could not fetch nf-core repositories! $gh_repos_url");
}
foreach($gh_repos as $repo){
    if(in_array($repo->name, $ignored_repos)){
        $repo_type = 'core_repos';
    } else {
        $repo_type = 'pipelines';
    }
    $results[$repo_type][$repo->name]['repo_metrics'][$updated] = array(
        'id' => $repo->id,
        'name' => $repo->name,
        'full_name' => $repo->full_name,
        'private' => $repo->private,
        'html_url' => $repo->html_url,
        'description' => $repo->description,
        'created_at' => $repo->created_at,
        'updated_at' => $repo->updated_at,
        'pushed_at' => $repo->pushed_at,
        'size' => $repo->size,
        'stargazers_count' => $repo->stargazers_count,
        'forks_count' => $repo->forks_count,
        'archived' => $repo->archived
    );
}

// Fetch new statistics for each repo
foreach(['pipelines', 'core_repos'] as $repo_type){
    foreach($results[$repo_type] as $repo_name => $repo_stats){
        // Views
        $gh_views_url = 'https://api.github.com/repos/nf-core/'.$repo_name.'/traffic/views';
        $gh_views = json_decode(file_get_contents($gh_views_url, false, $gh_api_opts));
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core repo views! $gh_views_url");
        }
        foreach($gh_views->views as $view){
            $results[$repo_type][$repo_name]['views_count'][$view->timestamp] = $view->count;
            $results[$repo_type][$repo_name]['views_uniques'][$view->timestamp] = $view->uniques;
        }
        // Clones
        $gh_clones_url = 'https://api.github.com/repos/nf-core/'.$repo_name.'/traffic/clones';
        $gh_clones = json_decode(file_get_contents($gh_clones_url, false, $gh_api_opts));
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core repo clones! $gh_clones_url");
        }
        foreach($gh_clones->clones as $clone){
            $results[$repo_type][$repo_name]['clones_count'][$clone->timestamp] = $clone->count;
            $results[$repo_type][$repo_name]['clones_uniques'][$clone->timestamp] = $clone->uniques;
        }
        // Contributors
        $gh_contributors_url = 'https://api.github.com/repos/nf-core/'.$repo_name.'/stats/contributors';
        $gh_contributors = json_decode(file_get_contents($gh_contributors_url, false, $gh_api_opts));
        // If the data hasn't been cached when you query a repository's statistics, you'll receive a 202 response;
        // a background job is also fired to start compiling these statistics.
        // Give the job a few moments to complete, and then submit the request again
        if(in_array("HTTP/1.1 202 Accepted", $http_response_header)){
            $contribs_try_again[$repo_name] = $gh_contributors_url;
        } else if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core repo contributors! $gh_contributors_url");
        }
        $results[$repo_type][$repo_name]['contributors'] = $gh_contributors;
        $results[$repo_type][$repo_name]['num_contributors'] = count($gh_contributors);

        // Recalculate totals
        foreach(['views_count', 'views_uniques', 'clones_count', 'clones_uniques'] as $ctype){
            $results[$repo_type][$repo_name][$ctype.'_total'] = 0;
            if(count($results[$repo_type][$repo_name][$ctype]) > 0){
                foreach($results[$repo_type][$repo_name][$ctype] as $stat){
                    $results[$repo_type][$repo_name][$ctype.'_total'] += $stat;
                }
            }
        }
    }
}

// Try contribs again now that we've let it fire
if(count($contribs_try_again) > 0){
    sleep(10);
    foreach($contribs_try_again as $repo_name => $gh_contributors_url){
        $gh_contributors = json_decode(file_get_contents($gh_contributors_url, false, $gh_api_opts));
        if(in_array("HTTP/1.1 202 Accepted", $http_response_header)){
            echo("Tried getting contributors after delay for $repo_name, but took too long.");
        } else if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            var_dump($http_response_header);
            die("Could not fetch nf-core repo contributors! $gh_contributors_url");
        }
        $results[$repo_type][$repo_name]['contributors'] = $gh_contributors;
        $results[$repo_type][$repo_name]['num_contributors'] = count($gh_contributors);
    }
}

// The data for commits per week is massive - remove it
foreach(['pipelines', 'core_repos'] as $repo_type){
    foreach($results[$repo_type] as $repo_name => $repo_stats){
        foreach($results[$repo_type][$repo_name]['contributors'] as $idx => $contributor){
            unset($results[$repo_type][$repo_name]['contributors'][$idx]->weeks);
        }
    }
}



//
//
// SLACK USERS
//
//

$slack_api_url = 'https://slack.com/api/team.billableInfo?token='.$config['slack_access_token'].'&pretty=1';
$slack_api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => [
            'User-Agent: PHP',
            'Content-Type: application/x-www-form-urlencoded'
        ]
    ]
]);
$slack_users = json_decode(file_get_contents($slack_api_url, false, $slack_api_opts));
if(!in_array("HTTP/1.1 200 OK", $http_response_header) || !isset($slack_users->ok) || !$slack_users->ok){
    var_dump($http_response_header);
    echo("Could not fetch slack user list!");
} else {
    $results['slack']['user_counts'][$updated] = [
        'total' => 0,
        'active' => 0,
        'inactive' => 0
    ];
    foreach($slack_users->billable_info as $uid => $user){
        $results['slack']['user_counts'][$updated]['total'] += 1;
        if($user->billing_active){
            $results['slack']['user_counts'][$updated]['active'] += 1;
        } else {
            $results['slack']['user_counts'][$updated]['inactive'] += 1;
        }
    }
}



//
//
// DONE - save results to JSON file
//
//

// Print results to a file
$results_json = json_encode($results, JSON_PRETTY_PRINT)."\n";
file_put_contents($results_fn, $results_json);

echo("update_stats done " . mktime()."\n");
