<?php
// ini_set('display_errors', 1);
// ini_set('display_startup_errors', 1);
// error_reporting(E_ALL);

// start the clocks
$start_time = microtime(TRUE);
require_once('../includes/functions.php');


if (!$curr_event) {
    $title = 'Sorry,';
    $subtitle = 'there is no event on-going. 😞';
    include('../includes/header.php');
} else {
    $title = $curr_event['title'];
    $subtitle = 'Measuring activity across the nf-core community during this event.';
    $import_chartjs = true;
    echo '<script>var theme = "' . $theme . '"</script>';
    include('../includes/header.php');

    $config = parse_ini_file("../config.ini");

    $conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);
    if ($conn === false) {
        echo "<script>console.log('No records matching the query " . mysqli_connect_error() . " were found.');</script>";
    }

    // Attempt select query execution
    $sql = "SELECT * FROM github_org_members";
    if ($result = mysqli_query($conn, $sql)) {
        if (mysqli_num_rows($result) > 0) {
            $gh_members_num = mysqli_num_rows($result);
            // Free result set
            mysqli_free_result($result);
        } else {
            echo "<script>console.log('No records matching the query " . $sql . " were found.');</script>";
        }
    } else {
        echo "<script>console.log(`ERROR: Was not able to execute " . json_encode($sql) . json_encode(mysqli_error($conn)) . "`);</script>";
    }

    // Attempt select query execution
    $sql = "SELECT * FROM slack_users";
    if ($result = mysqli_query($conn, $sql)) {
        if (mysqli_num_rows($result) > 0) {
            $slack_users_num = mysqli_num_rows($result);
            // Free result set
            mysqli_free_result($result);
        } else {
            echo "<script>console.log('No records matching the query " . $sql . " were found.');</script>";
        }
    } else {
        echo "<script>console.log(`ERROR: Was not able to execute " . json_encode($sql) . json_encode(mysqli_error($conn)) . "`);</script>";
    }

    function github_query($gh_query_url)
    {
        global $config;
        $gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

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
            if (strpos($http_response_header[0], "HTTP/1.1 200") === false) {
                var_dump($http_response_header);
                echo ("\nCould not fetch $gh_query_url");
                continue;
            }

            array_push($res, ...$tmp_results);

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
?>

    <h1 id="community"><a href="#community" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Community</h1>
    <p>The numbers below track our growth over the various channels that the nf-core community operates in.</p>
    <p class="text-info small">
        <i class="far fa-hand-point-right"></i>
        Click a number to see how the community has grown over time
    </p>
    <div class="card-group text-center stats_keynumbers">
        <div class="card bg-light">
            <div class="card-body">
                <p class="card-text display-4"><a href="#gh_orgmembers" class="text-body text-decoration-none stretched-link"><?php echo $gh_members_num; ?></a></p>
                <p class="card-text text-muted">GitHub organisation members</p>
            </div>
            <div class="bg-icon"><i class="fab fa-github"></i></div>
        </div>
        <div class="card bg-light">
            <div class="card-body">
                <p class="card-text display-4"><a href="#slack_users" class="text-body text-decoration-none stretched-link"><?php echo $slack_users_num; ?></a></p>
                <p class="card-text text-muted">Slack users</p>
            </div>
            <div class="bg-icon"><i class="fab fa-slack"></i></div>
        </div>
</div>
    <?php

} // close else-statement from the very beginning


// Stop the clocks!
$end_time = microtime(TRUE);
$time_taken = round($end_time - $start_time, 5);

$subfooter = '<p class="mb-0"><i class="far fa-clock"></i> Last updated: ' . date('d-m-Y', $stats_json->updated) . '. Page generated in ' . $time_taken . ' seconds.</p>';



include('../includes/footer.php'); ?>