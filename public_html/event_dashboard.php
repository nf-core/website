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

    require_once '../dbconfig.php';

    $conn = mysqli_connect($host, $username, $password, $dbname, $port);

    // Attempt select query execution
    $sql = "SELECT * FROM github_org_users";
    if ($result = mysqli_query($conn, $sql)) {
        if (mysqli_num_rows($result) > 0) {
            $members_num = mysqli_num_rows($result);
            // Free result set
            mysqli_free_result($result);
        } else {
            echo "<script>console.log('No records matching the query " . $sql . " were found.');</script>";
        }
    } else {
        echo "<script>console.log(`ERROR: Was not able to execute " . json_encode($sql) . json_encode(mysqli_error($conn)) . "`);</script>";
    }

?>

    <h1 id="community"><a href="#community" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Community</h1>
    <p>The numbers below track our growth over the various channels that the nf-core community operates in.</p>
    <p class="text-info small">
        <i class="far fa-hand-point-right"></i>
        Click a number to see how the community has grown over time
    </p>
    <div class="card bg-light">
        <div class="card-body">
            <p class="card-text display-4"><a href="#gh_orgmembers" class="text-body text-decoration-none stretched-link"><?php echo $members_num; ?></a></p>
            <p class="card-text text-muted">GitHub organisation members</p>
        </div>
        <div class="bg-icon"><i class="fab fa-github"></i></div>
    </div>
    

<?php

} // close else-statement from the very beginning


// Stop the clocks!
$end_time = microtime(TRUE);
$time_taken = round($end_time - $start_time, 5);

$subfooter = '<p class="mb-0"><i class="far fa-clock"></i> Last updated: ' . date('d-m-Y', $stats_json->updated) . '. Page generated in ' . $time_taken . ' seconds.</p>';



include('../includes/footer.php'); ?>