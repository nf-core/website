<?php
// Import entries from nfcore_stats.json and nfcore_issue_stats.json into the database

// IMPORTANT! the corresponding tables have to exist. Be sure to run update_stats.php at least once before this script.

echo "\nRunning import_stats_json - " . date('Y-m-d h:i:s') . "\n";
$config = parse_ini_file('config.ini');
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

if ($conn === false) {
    die('ERROR: Could not connect. ' . mysqli_connect_error());
}

// get all pipelines
$sql = 'SELECT * FROM nfcore_pipelines';
$pipelines = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipelines = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}

$stats_json_fn = dirname(__FILE__) . '/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn), true);

$issues_json_fn = dirname(__FILE__) . '/nfcore_issue_stats.json';
$issues_json = json_decode(file_get_contents($issues_json_fn), true);

// Prepare an insert statement
$sql =
    'INSERT INTO github_traffic_stats ( pipeline_id,views,views_uniques,clones,clones_uniques,timestamp) VALUES (?,?,?,?,?,?)';

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, 'iiiiis', $pipeline_id, $views, $views_uniques, $clones, $clones_uniques, $timestamp);

    foreach ($pipelines as $idx => $pipeline) {
        $gh_views = $stats_json['pipelines'][$pipeline['name']]['views_count'];
        if($gh_views){
            foreach ($gh_views as $timestamp_raw => $views_count) {
                $timestamp = date('Y-m-d H:i:s', strtotime($timestamp_raw));
                $check =
                    "SELECT * FROM github_traffic_stats WHERE pipeline_id = '" .
                    $pipeline['id'] .
                    "' AND timestamp = '" .
                    $timestamp . "'";
                $res = mysqli_query($conn, $check);
                if ($res->num_rows) {
                    echo 'Entry for pipeline ' .
                        $pipeline['bane'] .
                        " and timestamp '" .
                        $timestamp .
                        "' already exists. Skipping.\n";
                    continue;
                } else {
                    $pipeline_id = $pipeline['id'];
                    $views = $views_count;
                    $views_uniques = $stats_json['pipelines'][$pipeline['name']]['views_uniques'][$timestamp_raw];
                    $clones = $stats_json['pipelines'][$pipeline['name']]['clones_count'][$timestamp_raw];
                    $clones_uniques = $stats_json['pipelines'][$pipeline['name']]['clones_uniques'][$timestamp_raw];

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
