<?php
//
// ---------------------------
// GitHub shows traffic in the form of repo views and clones, however
// the data is only available for two weeks.
// We want it forever! So this script scrapes and saves the data.
// It is intended to be run routinely using a cronjob
//
// Manual usage: on command line, simply execute this script:
//   $ php update_stats.php

// Allow PHP fopen to work with remote links
ini_set('allow_url_fopen', 1);

echo "\nRunning update_stats - " . date('Y-m-d h:i:s') . "\n";
$config = parse_ini_file('config.ini');
$gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

if ($conn === false) {
    die('ERROR: Could not connect. ' . mysqli_connect_error());
}

// get all pipelines
$sql = 'SELECT * FROM nfcore_pipelines WHERE pipeline_type = "pipelines"';
$pipelines = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipelines = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}

function github_query($gh_query_url) {
    global $config;
    $gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

    // HTTP header to use on GitHub API GET requests
    $gh_api_opts = stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => ['User-Agent: PHP', "Authorization: Basic $gh_auth"],
        ],
    ]);

    $first_page = true;
    $next_page = false;
    $res = [];
    while ($first_page || $next_page) {
        // reset loop vars
        $first_page = false;
        // Get GitHub API results
        if ($next_page) {
            $gh_query_url = $next_page;
        }
        $tmp_results = json_decode(file_get_contents($gh_query_url, false, $gh_api_opts), true);
        // If the data hasn't been cached when you query a repository's statistics, you'll receive a 202 response;
        // a background job is also fired to start compiling these statistics.
        // Give the job a few moments to complete, and then submit the request again
        if (strpos($http_response_header[0], 'HTTP/1.1 202') !== false) {
            echo "\nWaiting for GitHub API to return results for $gh_query_url";
            sleep(5);
            $first_page = true;
            continue;
        }
        if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
            var_dump($http_response_header);
            echo "\nCould not fetch $gh_query_url";
            continue;
        }
        if (substr($gh_query_url, 0, 29) == 'https://api.github.com/repos/') {
            $res = $tmp_results;
        } else {
            array_push($res, ...$tmp_results);
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
    return $res;
}

if (!mysqli_query($conn, $sql)) {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}
// create github_pipeline_contrib_stats table with foreign keys to pipelines
$sql = "CREATE TABLE IF NOT EXISTS github_pipeline_contrib_stats (
        id             INT             AUTO_INCREMENT PRIMARY KEY,
        pipeline_id    INT             NOT NULL,
        author         VARCHAR(255)    NOT NULL,
        avatar_url     VARCHAR(255)    NOT NULL,
        week_date      datetime        NOT NULL,
        week_additions INT             NOT NULL,
        week_deletions INT             NOT NULL,
        week_commits   INT             NOT NULL,
        FOREIGN KEY (pipeline_id)   REFERENCES nfcore_pipelines(id)
        )";
if (mysqli_query($conn, $sql)) {
    echo "`github_pipeline_contrib_stats` table created successfully.\n";
} else {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}

// Prepare an insert statement
$sql =
    'INSERT INTO github_pipeline_contrib_stats (pipeline_id, author, avatar_url, week_date, week_additions, week_deletions, week_commits) VALUES (?, ?, ?, ?, ?, ?, ?)';

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param(
        $stmt,
        'isssiii',
        $pipeline_id,
        $author,
        $avatar_url,
        $week_date,
        $week_additions,
        $week_deletions,
        $week_commits,
    );
    // add modules repo to pipelines
    foreach ($pipelines as $idx => $pipeline) {
        echo "\nUpdating stats for " . $pipeline['name'];
        // get contributors
        $gh_contributors = github_query(
            'https://api.github.com/repos/nf-core/' . $pipeline['name'] . '/stats/contributors',
        );
        foreach ($gh_contributors as $contributor) {
            $pipeline_id = $pipeline['id'];
            $author = $contributor['author']['login'];
            $avatar_url = $contributor['author']['avatar_url'];
            foreach ($contributor['weeks'] as $week) {
                $week_date = date('Y-m-d', $week['w']);
                //check if entry for this pipeline_id and week_date already exists and skip if so
                $check =
                    "SELECT * FROM github_pipeline_contrib_stats WHERE pipeline_id = '" .
                    $pipeline_id .
                    "' AND author = '" .
                    $author .
                    "' AND week_date = '" .
                    $week_date .
                    "'";
                $res = mysqli_query($conn, $check);
                if ($res->num_rows) {
                    continue;
                } else {
                    $week_additions = $week['a'];
                    $week_deletions = $week['d'];
                    $week_commits = $week['c'];
                    if (!mysqli_stmt_execute($stmt)) {
                        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
                    }
                }
            }
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

if (!mysqli_query($conn, $sql)) {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}
// create github_traffic_stats table with foreign keys to pipelines
$sql = "CREATE TABLE IF NOT EXISTS github_traffic_stats (
            id             INT             AUTO_INCREMENT PRIMARY KEY,
            pipeline_id    INT             NOT NULL,
            views          INT             DEFAULT NULL,
            views_uniques  INT             DEFAULT NULL,
            clones         INT             DEFAULT NULL,
            clones_uniques INT             DEFAULT NULL,
            timestamp      datetime        NOT NULL,
            FOREIGN KEY (pipeline_id)   REFERENCES nfcore_pipelines(id)
            )";
if (mysqli_query($conn, $sql)) {
    echo "`github_traffic_stats` table created successfully.\n";
} else {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}

// Prepare an insert statement
$sql =
    'INSERT INTO github_traffic_stats ( pipeline_id,views,views_uniques,clones,clones_uniques,timestamp) VALUES (?,?,?,?,?,?)';

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, 'iiiiis', $pipeline_id, $views, $views_uniques, $clones, $clones_uniques, $timestamp);

    foreach ($pipelines as $idx => $pipeline) {
        echo "\nUpdating traffic stats for " . $pipeline['name'];
        $gh_views = github_query('https://api.github.com/repos/nf-core/' . $pipeline['name'] . '/traffic/views');
        $gh_clones = github_query('https://api.github.com/repos/nf-core/' . $pipeline['name'] . '/traffic/clones');
        foreach ($gh_views['views'] as $gh_view) {
            $timestamp = date('Y-m-d H:i:s', strtotime($gh_view['timestamp']));
            $check =
                "SELECT * FROM github_traffic_stats WHERE pipeline_id = '" .
                $pipeline['id'] .
                "' AND timestamp = '" .
                $timestamp .
                "'";
            $res = mysqli_query($conn, $check);
            if ($res->num_rows) {
                echo "\n Entry already exists for pipeline_id " .
                    $pipeline['id'] .
                    ' and timestamp ' .
                    $timestamp .
                    ' db_timestamp ';
                continue;
            } else {
                // get gh_clones where timestamp matches
                foreach ($gh_clones['clones'] as $gh_clone) {
                    if ($gh_clone['timestamp'] == $gh_view['timestamp']) {
                        $gh_clones_count = $gh_clone['count'];
                        $gh_clones_unique = $gh_clone['uniques'];
                        break;
                    }
                }
                $pipeline_id = $pipeline['id'];
                $views = $gh_view['count'];
                $views_uniques = $gh_view['uniques'];
                $clones = $gh_clones_count;
                $clones_uniques = $gh_clones_unique;

                if (!mysqli_stmt_execute($stmt)) {
                    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
                }
            }
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}
mysqli_close($conn);
echo "\n Finished updating the database - " . date('Y-m-d h:i:s') . "\n";

###########################################################################################################
#                                                                                                         #
#                                              🕱 OLD CODE 🕱                                              #
#                                                                                                         #
###########################################################################################################

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
//   $ php update_stats.php

// Allow PHP fopen to work with remote links
ini_set('allow_url_fopen', 1);

// Use same updated time for everything
$updated = time();

// Final filename to write JSON to
$results_fn = dirname(__FILE__) . '/nfcore_stats.json';
$contribs_fn_root = dirname(__FILE__) . '/contributor_stats/';

echo "\nRunning update_stats - " . date('Y-m-d h:i:s') . "\n";

// Initialise the results array with the current time and placeholders
$results = [
    'updated' => $updated,
    'pipelines' => [],
    'core_repos' => [],
    'slack' => [],
    'gh_org_members' => [],
    'gh_contributors' => [],
    'gh_commits' => [],
    'gh_additions' => [],
    'gh_deletions' => [],
];

// Load a copy of the existing JSON file, if it exists
if (file_exists($results_fn)) {
    $results = json_decode(file_get_contents($results_fn), true);
}
$results['updated'] = $updated;

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
$pipelines_json = json_decode(file_get_contents(dirname(__FILE__) . '/public_html/pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
$contribs_try_again = [];

// Build array of repos to query
$pipelines_json_names = [];
foreach ($pipelines as $wf) {
    $pipelines_json_names[] = $wf->name;
    if (!isset($results['pipelines'][$wf->name])) {
        $results['pipelines'][$wf->name] = [];
    }
    $results['pipelines'][$wf->name]['num_releases'] = count($wf->releases);
}
$ignored_repos = parse_ini_file('ignored_repos.ini')['repos'];
foreach ($ignored_repos as $name) {
    if (!isset($results['core_repos'][$name])) {
        $results['core_repos'][$name] = [];
    }
}

// Delete cached pipelines stats for pipelines that have been deleted
foreach (array_keys($results['pipelines']) as $wfname) {
    if (!in_array($wfname, $pipelines_json_names)) {
        echo "Removing $wfname from the cached results as it appears to have been deleted.\n";
        unset($results['pipelines'][$wfname]);
    }
}

// Get snapshot of key metrics for all repos

// Get the current number of organisation members
// Returns 30 results per page!
$gh_members_url = 'https://api.github.com/orgs/nf-core/members';
$results['gh_org_members'][$updated] = 0;
$first_page = true;
$next_page = false;
while ($first_page || $next_page) {
    // reset loop vars
    $first_page = false;
    // Get GitHub API results
    if ($next_page) {
        $gh_members_url = $next_page;
    }
    $gh_members = json_decode(file_get_contents($gh_members_url, false, $gh_api_opts));
    if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
        var_dump($http_response_header);
        echo "Could not fetch nf-core members! $gh_members_url";
        continue;
    }
    $results['gh_org_members'][$updated] += count($gh_members);
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

// Fetch all repositories at nf-core
// $gh_repos_url = 'https://api.github.com/orgs/nf-core/repos?per_page=100';
// $gh_repos = json_decode(file_get_contents($gh_repos_url, false, $gh_api_opts));
// if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//     var_dump($http_response_header);
//     die("Could not fetch nf-core repositories! $gh_repos_url");
// }
// foreach ($gh_repos as $repo) {
//     if (in_array($repo->name, $ignored_repos)) {
//         $repo_type = 'core_repos';
//     } else {
//         $repo_type = 'pipelines';
//     }
//     $results[$repo_type][$repo->name]['repo_metrics'][$updated] = [
//         'id' => $repo->id,
//         'name' => $repo->name,
//         'full_name' => $repo->full_name,
//         'private' => $repo->private,
//         'html_url' => $repo->html_url,
//         'description' => $repo->description,
//         'created_at' => $repo->created_at,
//         'updated_at' => $repo->updated_at,
//         'pushed_at' => $repo->pushed_at,
//         'size' => $repo->size,
//         'stargazers_count' => $repo->stargazers_count,
//         'forks_count' => $repo->forks_count,
//         'archived' => $repo->archived,
//     ];
//     // Annoyingly, two values are only available if we query for just this repo
//     $gh_repo_url = 'https://api.github.com/repos/nf-core/' . $repo->name;
//     $gh_repo = json_decode(file_get_contents($gh_repo_url, false, $gh_api_opts));
//     if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//         var_dump($http_response_header);
//         echo "Could not fetch nf-core repo! $gh_repo_url";
//         continue;
//     }
//     $results[$repo_type][$repo->name]['repo_metrics'][$updated]['network_forks_count'] = $gh_repo->network_count;
//     $results[$repo_type][$repo->name]['repo_metrics'][$updated]['subscribers_count'] = $gh_repo->subscribers_count;
// }

// // Fetch new statistics for each repo
// foreach (['pipelines', 'core_repos'] as $repo_type) {
//     foreach ($results[$repo_type] as $repo_name => $repo_stats) {
//         // Views
//         $gh_views_url = 'https://api.github.com/repos/nf-core/' . $repo_name . '/traffic/views';
//         $gh_views = json_decode(file_get_contents($gh_views_url, false, $gh_api_opts));
//         if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//             // Pipelines are removed from the cache earlier as we know their names
//             if ($repo_type == 'core_repos' && strpos($http_response_header[0], 'HTTP/1.1 404') !== false) {
//                 echo 'Removing ' . $repo_name . " from the cached results as it appears to have been deleted.\n";
//                 unset($results['core_repos'][$repo_name]);
//             } else {
//                 echo "--------   Could not fetch nf-core repo views! $gh_views_url\n";
//                 var_dump($http_response_header);
//                 echo "\n--------   End of header for $gh_views_url\n\n\n";
//             }
//             continue;
//         }
//         foreach ($gh_views->views as $view) {
//             $results[$repo_type][$repo_name]['views_count'][$view->timestamp] = $view->count;
//             $results[$repo_type][$repo_name]['views_uniques'][$view->timestamp] = $view->uniques;
//         }
//         // Clones
//         $gh_clones_url = 'https://api.github.com/repos/nf-core/' . $repo_name . '/traffic/clones';
//         $gh_clones = json_decode(file_get_contents($gh_clones_url, false, $gh_api_opts));
//         if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//             var_dump($http_response_header);
//             echo "Could not fetch nf-core repo clones! $gh_clones_url";
//             continue;
//         }
//         foreach ($gh_clones->clones as $clone) {
//             $results[$repo_type][$repo_name]['clones_count'][$clone->timestamp] = $clone->count;
//             $results[$repo_type][$repo_name]['clones_uniques'][$clone->timestamp] = $clone->uniques;
//         }
//         // Contributors
//         $gh_contributors_url = 'https://api.github.com/repos/nf-core/' . $repo_name . '/stats/contributors';
//         $gh_contributors_raw = file_get_contents($gh_contributors_url, false, $gh_api_opts);
//         file_put_contents($contribs_fn_root . $repo_name . '.json', $gh_contributors_raw);
//         $gh_contributors = json_decode($gh_contributors_raw);
//         // If the data hasn't been cached when you query a repository's statistics, you'll receive a 202 response;
//         // a background job is also fired to start compiling these statistics.
//         // Give the job a few moments to complete, and then submit the request again
//         if (strpos($http_response_header[0], 'HTTP/1.1 202') !== false) {
//             $contribs_try_again[$repo_name] = [
//                 'repo_type' => $repo_type,
//                 'gh_contributors_url' => $gh_contributors_url,
//             ];
//         } elseif (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//             var_dump($http_response_header);
//             echo "Could not fetch nf-core repo contributors! $gh_contributors_url";
//             continue;
//         }
//         $results[$repo_type][$repo_name]['contributors'] = $gh_contributors;
//         $results[$repo_type][$repo_name]['num_contributors'] = count($gh_contributors);

//         // Commits
//         $results[$repo_type][$repo_name]['commits'] = 0;
//         foreach ($gh_contributors as $contributor) {
//             $results[$repo_type][$repo_name]['commits'] += $contributor->total;
//         }

//         // Recalculate totals
//         foreach (['views_count', 'views_uniques', 'clones_count', 'clones_uniques'] as $ctype) {
//             $results[$repo_type][$repo_name][$ctype . '_total'] = 0;
//             if (
//                 isset($results[$repo_type][$repo_name][$ctype]) &&
//                 count($results[$repo_type][$repo_name][$ctype]) > 0
//             ) {
//                 foreach ($results[$repo_type][$repo_name][$ctype] as $stat) {
//                     $results[$repo_type][$repo_name][$ctype . '_total'] += $stat;
//                 }
//             }
//         }
//     }
// }

// // Try contribs again now that we've let it fire
// if (count($contribs_try_again) > 0) {
//     sleep(5);
//     foreach ($contribs_try_again as $repo_name => $details) {
//         extract($details); // $repo_type, $gh_contributors_raw
//         $gh_contributors_raw = file_get_contents($gh_contributors_url, false, $gh_api_opts);
//         file_put_contents($contribs_fn_root . $repo_name . '.json', $gh_contributors_raw);
//         $gh_contributors = json_decode($gh_contributors_raw);
//         if (strpos($http_response_header[0], 'HTTP/1.1 202') !== false) {
//             echo "Tried getting contributors after delay for $repo_name, but took too long.";
//             continue;
//         } elseif (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
//             var_dump($http_response_header);
//             echo "Could not fetch nf-core repo contributors! $gh_contributors_url";
//             continue;
//         }
//         $results[$repo_type][$repo_name]['contributors'] = $gh_contributors;
//         $results[$repo_type][$repo_name]['num_contributors'] = count($gh_contributors);
//     }
// }

// foreach (['pipelines', 'core_repos'] as $repo_type) {
//     foreach ($results[$repo_type] as $repo_name => $repo_stats) {
//         foreach ($results[$repo_type][$repo_name]['contributors'] as $idx => $contributor) {
//             // Count how many total contributors and contributions we have per week
//             foreach ($contributor->weeks as $w) {
//                 // Skip zeros (anything before 2010)
//                 if ($w->w < 1262304000) {
//                     continue;
//                 }
//                 // Find earliest contribution per author
//                 if (!isset($results['gh_contributors'][$contributor->author->login])) {
//                     $results['gh_contributors'][$contributor->author->login] = $w->w;
//                 }
//                 $results['gh_contributors'][$contributor->login] = min(
//                     $w->w,
//                     $results['gh_contributors'][$contributor->login],
//                 );
//                 // Sum total contributions for everyone
//                 if (!isset($results['gh_commits'][$w->w])) {
//                     $results['gh_commits'][$w->w] = 0;
//                 }
//                 if (!isset($results['gh_additions'][$w->w])) {
//                     $results['gh_additions'][$w->w] = 0;
//                 }
//                 if (!isset($results['gh_deletions'][$w->w])) {
//                     $results['gh_deletions'][$w->w] = 0;
//                 }
//                 $results['gh_commits'][$w->w] += $w->c;
//                 $results['gh_additions'][$w->w] += $w->a;
//                 $results['gh_deletions'][$w->w] += $w->d;
//             }

//             // The data for commits per week is massive - remove it
//             unset($results[$repo_type][$repo_name]['contributors'][$idx]->weeks);
//         }
//     }
// }

//
//
// SLACK USERS
//
//
echo 'update_stats - Slack updates - ' . date('Y-m-d h:i:s') . "\n";

$slack_api_url = 'https://slack.com/api/team.billableInfo?token=' . $config['slack_access_token'] . '&pretty=1';
$slack_api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => ['User-Agent: PHP', 'Content-Type: application/x-www-form-urlencoded'],
    ],
]);
$slack_users = json_decode(file_get_contents($slack_api_url, false, $slack_api_opts));
if (strpos($http_response_header[0], 'HTTP/1.0 200') === false || !isset($slack_users->ok) || !$slack_users->ok) {
    var_dump($http_response_header);
    echo 'Could not fetch slack user list!';
} else {
    $results['slack']['user_counts'][$updated] = [
        'total' => 0,
        'active' => 0,
        'inactive' => 0,
    ];
    foreach ($slack_users->billable_info as $uid => $user) {
        $results['slack']['user_counts'][$updated]['total'] += 1;
        if ($user->billing_active) {
            $results['slack']['user_counts'][$updated]['active'] += 1;
        } else {
            $results['slack']['user_counts'][$updated]['inactive'] += 1;
        }
    }
}

// //
// //
// // Twitter - get number of followers
// //
// //
// echo 'update_stats - twitter updates - ' . date('Y-m-d h:i:s') . "\n";
// require 'vendor/autoload.php';
// use Abraham\TwitterOAuth\TwitterOAuth;
// // Connect to twitter
// $connection = new TwitterOAuth(
//     $config['twitter_key'],
//     $config['twitter_secret'],
//     $config['twitter_access_token'],
//     $config['twitter_access_token_secret'],
// );
// $twitter_stats = $connection->get('users/show', ['screen_name' => 'nf_core']);
// if (isset($twitter_stats->followers_count)) {
//     $results['twitter']['followers_count'][$updated] = $twitter_stats->followers_count;
// }

//
//
// DONE - save results to JSON file
//
//

// Print results to a file
echo "update_stats - Saving to $results_fn - " . date('Y-m-d h:i:s') . "\n";
$results_json = json_encode($results, JSON_PRETTY_PRINT) . "\n";
file_put_contents($results_fn, $results_json);

echo 'update_stats done ' . date('Y-m-d h:i:s') . "\n\n";
