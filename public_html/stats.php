<?php
// ini_set('display_errors', 1);
// ini_set('display_startup_errors', 1);
// error_reporting(E_ALL);
ob_start();
// Please someone rewrite this page..
set_time_limit(240);

// start the clocks
$start_time = microtime(true);

$title = 'nf-core in numbers';
$subtitle = 'Measuring activity across the nf-core community.';
$import_chartjs = true;
include '../includes/slim_header.php';

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

$stats_json_fn = dirname(dirname(__FILE__)) . '/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn));

$issues_json_fn = dirname(dirname(__FILE__)) . '/nfcore_issue_stats.json';
$issues_json = json_decode(file_get_contents($issues_json_fn), true);
$issues_json_latest = false;
foreach ($issues_json['stats'] as $key => $arr) {
    if (!is_numeric($key)) {
        continue;
    }
    if (!$issues_json_latest) {
        $issues_json_latest = $key;
    }
    $issues_json_latest = max($issues_json_latest, $key);
}

// Convenience variables
$latest_slack_update = max(array_keys(get_object_vars($stats_json->slack->user_counts)));
$slack_users = $stats_json->slack->user_counts->{$latest_slack_update};
$twitter_datekeys = array_keys(get_object_vars($stats_json->twitter->followers_count));
$twitter_users = $stats_json->twitter->followers_count->{max($twitter_datekeys)};

$config = parse_ini_file('../config.ini');
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);
// get stats for pipelines
$sql = 'SELECT timestamp,
                SUM(views) AS sum_total_views,
                SUM(views_uniques) AS sum_total_views_uniques,
                SUM(clones) AS sum_total_clones,
                SUM(clones_uniques) AS sum_total_clones_uniques,
                MIN(timestamp) AS min_timestamp FROM github_traffic_stats
            INNER JOIN nfcore_pipelines ON github_traffic_stats.pipeline_id = nfcore_pipelines.id GROUP BY timestamp ORDER BY timestamp ASC';
$pipeline_metrics = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipeline_metrics = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}
$view_counts = [];
$clone_counts = [];
# get views per day
$sql = 'SELECT timestamp, SUM(views) AS sum_total_views, SUM(views_uniques) AS sum_total_views_unique FROM github_traffic_stats
            INNER JOIN nfcore_pipelines ON github_traffic_stats.pipeline_id = nfcore_pipelines.id GROUP BY timestamp ORDER BY timestamp ASC';
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $view_counts = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}
$sql = 'SELECT timestamp, SUM(clones) AS sum_total_clones, SUM(clones_uniques) AS sum_total_clones_unique FROM github_traffic_stats
            INNER JOIN nfcore_pipelines ON github_traffic_stats.pipeline_id = nfcore_pipelines.id GROUP BY timestamp ORDER BY timestamp ASC';
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $clone_counts = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}
// get all repos
$sql = 'SELECT * FROM nfcore_pipelines';
$gh_repos = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $gh_repos = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}
// get unique contributor names from github_pipeline_contrib_stats table from the mysql database
$sql = 'SELECT DISTINCT author, avatar_url, SUM(week_commits) AS total_sum_commits
FROM github_pipeline_contrib_stats GROUP BY author,avatar_url  ORDER BY total_sum_commits DESC';
$gh_contributors_db = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $gh_contributors_db = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
} else {
    echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
}

// echo '<pre>' . print_r($gh_contributors_db, true) . '</pre>';
// get maximum total_sum_commits from $gh_contributors_db
$max_total_sum_commits = max(array_column($gh_contributors_db, 'total_sum_commits'));
$total_sum_commits = array_sum(array_column($gh_contributors_db, 'total_sum_commits'));
$metrics_since = min(array_column($pipeline_metrics, 'min_timestamp'));
$total_sum_clones = array_sum(array_column($pipeline_metrics, 'sum_total_clones'));
$total_sum_views = array_sum(array_column($pipeline_metrics, 'sum_total_views'));
$total_sum_clones_uniques = array_sum(array_column($pipeline_metrics, 'sum_total_clones_uniques'));
$total_sum_views_uniques = array_sum(array_column($pipeline_metrics, 'sum_total_views_uniques'));
# echo '<pre>'.print_r($stats, true).'</pre>';
// print_r($gh_contributors_db);
// Get unique contributors - commits and issues
$gh_contributors = [];
$gh_contributor_commits = [];
foreach( $gh_contributors_db as $contributor ) {
    $gh_contributor_commits[] = $contributor['author'];
    $gh_contributors[$contributor['author']] = $contributor['total_sum_commits'];
}
$gh_contributor_issues = array_keys($issues_json['authors']);
foreach ($issues_json['authors'] as $author => $info) {
        $gh_contributors[$author] = $info['first_contribution'];
}

?>

<h1>Introduction</h1>
<p>On this page you can see the beating heart of nf-core - the size of our community and the output of our work.</p>
<ul>
    <li><a href="#community">Community</a>
        <ul>
            <li><a href="#slack">Slack</a></li>
            <li><a href="#gh_orgmembers">GitHub organisation members</a></li>
            <li><a href="#gh_contribs">GitHub contributors</a></li>
            <li><a href="#twitter">Twitter followers</a></li>
        </ul>
    </li>
    <li><a href="#code">Code</a>
        <ul>
            <li><a href="#repo_traffic">Repository traffic</a></li>
            <li><a href="#github_prs">Pull Requests</a>, <a href="#github_pr_response_time">response times</a></li>
            <li><a href="#github_issues">Issues</a>, <a href="#github_issue_response_time">response times</a></li>
            <li><a href="#contributor_leaderboard">Contributor Leaderboard</a></li>
            <li><a href="#pipeline_numbers">Pipeline numbers</a></li>
            <li><a href="#pipelines">Pipelines</a></li>
            <li><a href="#core_repos">Core repositories</a></li>
        </ul>
    </li>
</ul>

<section id="caveats">
    <div class="card mb-3 mt-2">
        <div class="card-header">
            <a href="#caveats_<?php echo $repo_type; ?>" data-bs-toggle="collapse" data-bs-target="#caveats_<?php echo $repo_type; ?>" class="text-muted">
                <u>Click to expand:</u> How these numbers are collected and what caveats should be considered
            </a>
        </div>
        <div id="caveats_<?php echo $repo_type; ?>" class="collapse">
            <div class="card-body small">
                <p>Please bear in mind the following points when looking over these numbers:</p>
                <ul>
                    <li>Many pipelines are worked on long before they are forked to nf-core. The age, stars and other metrics of the original parent repository are not shown.</li>
                    <li>Metrics are for the default (<code>master</code>) branch only</li>
                    <li>Commits and contributors are only counted if associated with a GitHub account</li>
                    <li><code>nextflow pull</code> and <code>nextflow run</code> uses git to clone a remote repo the first time it runs, so the clones count gives some idea of usage. However:
                        <ul>
                            <li><em>Unique cloners</em> is based on IP address, so will under-represent institutional users sharing a single external IP address</li>
                            <li><em>Unique cloners</em> is based on IP address, so will over-represent cloud users using multiple IP addresses</li>
                            <li>Traditional HPC centres may share workflow installations, so only have one clone for many users &nbsp;/&nbsp; pipeline runs</li>
                            <li>Cloud users will typically spin up a new instance and clone the workflow every time that they run a pipeline.</li>
                        </ul>
                    </li>
                    <li>Clone counts and repositoriy views are only available for two weeks - longer term data collection for nf-core repos started in July 2019. This is when we started counting the totals.</li>
                    <li>Metrics are fetched once per day (last checked <?php echo date(
                        'Y-m-d',
                        $stats_json->updated,
                    ); ?>).</li>
                </ul>
            </div>
        </div>
    </div>
</section>

<?php echo _h2('Community'); ?>
<p>The numbers below track our growth over the various channels that the nf-core community operates in.</p>
<p class="text-info small">
    <i class="far fa-hand-point-right"></i>
    Click a number to see how the community has grown over time
</p>

<div class="card-group text-center stats_keynumbers">
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><a href="#slack" class="text-body text-decoration-none stretched-link"><?php echo $slack_users->total; ?></a></p>
            <p class="card-text text-muted">Slack users</p>
        </div>
        <div class="bg-icon"><i class="fab fa-slack"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><a href="#gh_orgmembers" class="text-body text-decoration-none stretched-link"><?php echo $stats_json
                ->gh_org_members->{$stats_json->updated}; ?></a></p>
            <p class="card-text text-muted">GitHub organisation members</p>
        </div>
        <div class="bg-icon"><i class="fab fa-github"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body" data-bs-toggle="tooltip" title="<?php echo count(
            $gh_contributor_commits,
        ); ?> have committed code, <?php echo count($gh_contributor_issues); ?> have written issues">
            <p class="card-text display-4"><a href="#gh_contribs" class="text-body text-decoration-none stretched-link"><?php echo count(
                $gh_contributors,
            ); ?></a></p>
            <p class="card-text text-muted">GitHub contributors</p>
        </div>
        <div class="bg-icon"><i class="fas fa-code-branch"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><a href="#twitter" class="text-body text-decoration-none stretched-link"><?php echo $twitter_users; ?></a></p>
            <p class="card-text text-muted">Twitter followers</p>
        </div>
        <div class="bg-icon"><i class="fab fa-twitter"></i></div>
    </div>
</div>

<div class="row">
    <div class="col-lg-6">
        <?php echo _h2('Slack'); ?>
        <p>Slack is a real-time messaging tool, with discussion split into channels and groups.
            We use it to provide help to people running nf-core pipelines, as well as discussing development ideas.
            You can join the nf-core slack by getting an invite <a href="https://nf-co.re/join/slack">here</a>.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="slack_users_plot" height="200"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> Slack considers users to be inactive when they haven't used slack for the previous 14 days.</p>
                <p class="card-text small text-muted mb-1"><i class="fas fa-exclamation-triangle"></i> Data from before 2019-07-24 fudged by reverse-engineering billing details on the slack admin pages.</p>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="slack" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="slack" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>

        <h2 class="mt-0" id="twitter">Twitter followers<a href="#twitter" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>We use our <a href="https://twitter.com/nf_core">@nf_core</a> twitter account to send automated tweets about new pipeline releases and other updates relevant to the community.
            Follower counts give some indication to the level of interest in the nf-core project.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="twitter_followers_plot" height="150"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-exclamation-triangle"></i> Data from before 2019-06-26 fudged by reverse-engineering a tiny sparkline plot on the twitter analytics website.</p>
                <p class="card-text small text-muted"><a href="#" data-bs-target="twitter" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="twitter" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
            </div>
        </div>
    </div>

    <div class="col-lg-6">

        <h2 class="mt-0" id="gh_orgmembers">GitHub organisation members<a href="#gh_orgmembers" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>We use GitHub to manage all of the code written for nf-core.
            It's a fantastic platform and provides a huge number of tools.
            We have a GitHub organisation called <a href="https://github.com/nf-core/">nf-core</a> which anyone can join:
            drop us a note <a href="https://github.com/nf-core/nf-co.re/issues/3">here</a> or anywhere and we'll send you an invite.
        </p>
        <p>It's not required to be a member of the nf-core GitHub organisation to contribute.
            However, members get the nf-core logo listed on their profile page and full write-access to all nf-core repositories.
        </p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="gh_orgmembers_plot" height="150"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-exclamation-triangle"></i> By default, organisation membership is private. This is why you'll see a lower number if you visit the <a href="https://github.com/nf-core/">nf-core organisation page</a> and are not a member.
                <p class="card-text small text-muted"><a href="#" data-bs-target="gh_orgmembers" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="gh_orgmembers" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
            </div>
        </div>

        <h2 class="mt-0" id="gh_contribs">GitHub contributors<a href="#gh_contribs" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>Anybody can fork nf-core repositories and open a pull-request.
            Here we count how many different people have contributed at least one commit to an nf-core repository, or created or commented on an issue or pull-request.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="gh_contribs_plot" height="180"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> Plot truncated to start of 2018 (some pipelines moved to nf-core so have older contributions).</p>
                <p class="card-text small text-muted"><a href="#" data-bs-target="gh_contribs" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="gh_contribs" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
            </div>
        </div>

    </div>
</div>


<h1 id="code">Code stats<a href="#code" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h1>
<p>Whilst we always prefer quality over quantity, these numbers reflect the work output from the nf-core community.</p>

<div class="card-group text-center stats_keynumbers">
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><?php echo count($pipelines) +
                count($pipelines); ?></p>
            <p class="card-text text-muted">Repositories</p>
        </div>
        <div class="bg-icon"><i class="far fa-folder"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><a href="#github_prs" class="text-body text-decoration-none stretched-link"><?php echo round_nicely(
                $issues_json['stats'][$issues_json_latest]['prs']['count'],
            ); ?></a></p>
            <p class="card-text text-muted">Pull Requests</p>
        </div>
        <div class="bg-icon"><i class="fas fa-code-branch fa-flip-vertical"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><?php echo round_nicely(
                $total_sum_commits,
            ); ?></p>
            <p class="card-text text-muted">Commits</p>
        </div>
        <div class="bg-icon"><i class="far fa-file-code"></i></div>
    </div>
    <div class="card bg-body">
        <div class="card-body">
            <p class="card-text display-4"><a href="#github_issues" class="text-body text-decoration-none stretched-link"><?php echo round_nicely(
                $issues_json['stats'][$issues_json_latest]['issues']['count'],
            ); ?></a></p>
            <p class="card-text text-muted">Issues</p>
        </div>
        <div class="bg-icon"><i class="fas fa-exclamation-circle"></i></div>
    </div>
</div>

<h2 class="mt-0" id="repo_traffic">Repository traffic<a href="#repo_traffic" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
<p>Every time a nextflow user pulls an nf-core pipeline, the repository is cloned. Here we can track how much that happens across all nf-core repositories.
    Please note that these numbers come with some caveats <a href="#caveats">[ see more ]</a>.</p>
<p>Additionally, GitHub tracks how many times people view repository web pages on github.com.</p>

<div class="card mt-4">
    <div class="card-header">
        <span class="float-end small text-muted">
            <a href="#" data-bs-target="repo_clones" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
            &nbsp;/&nbsp; <a href="#" data-bs-target="repo_clones" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
        </span>
        Git clones: All nf-core repositories
    </div>
    <div class="card-body">
        <canvas id="repo_clones_plot" height="80"></canvas>
    </div>
    <div class="card-footer text-muted text-center small">
        <div class="row">
            <div class="col-6 border-end border-secondary">
                <span class="text-body lead"><?php echo round_nicely(
                   $total_sum_clones,
                ); ?></span>
                <br>Clones since <?php
                    echo date('F Y', strtotime($metrics_since));
                ; ?>
            </div>
            <div class="col-6" data-bs-toggle="tooltip" title="Note: Unique per repository. Will double-count the same person cloning two different repositories.">
                <span class="text-body lead"><?php echo round_nicely(
                    $total_sum_clones_uniques,
                ); ?></span>
                <br>Unique cloners since <?php echo date('F Y', strtotime($metrics_since)); ?>
            </div>
        </div>
    </div>
</div>

<div class="card mt-4">
    <div class="card-header">
        <span class="float-end small text-muted">
            <a href="#" data-bs-target="repo_views" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
            &nbsp;/&nbsp; <a href="#" data-bs-target="repo_views" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
        </span>
        Visitors: All nf-core repositories
    </div>
    <div class="card-body">
        <canvas id="repo_views_plot" height="80"></canvas>
    </div>
    <div class="card-footer text-muted text-center small">
        <div class="row align-items-center">
            <div class="col-6 border-end border-secondary">
                <span class="text-body lead"><?php echo round_nicely(
                    $total_sum_views,
                ); ?></span>
                <br>Views since <?php echo date('F Y', strtotime($metrics_since)); ?>
            </div>
            <div class="col-6" data-bs-toggle="tooltip" title="Note: Unique per repository. Will double-count the same person viewing two different repositories.">
                <span class="text-body lead"><?php echo round_nicely(
                    $total_sum_views_uniques,
                ); ?></span>
                <br>Unique visitors since <?php echo date('F Y', strtotime($metrics_since)); ?>
            </div>
        </div>
    </div>
</div>

<div class="row">
    <div class="col-lg-6">

        <h2 class="mt-0" id="github_prs">Pull Requests<a href="#github_prs" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>When people contribute code to a nf-core repository, we conduct a "Pull request" - other members of the nf-core community review the proposed code and make suggestions, before merging into the main repository.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="github_prs_plot" height="200"></canvas>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="github_prs" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="github_prs" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>

        <h2 class="mt-0" id="github_pr_response_time">Pull Request response times<a href="#github_pr_response_time" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>Pull-requests are reviewed by the nf-core community - they can contain discussion on the code and can be merged and closed.
            We aim to be prompt with reviews and merging. Note that some PRs can be a simple type and so very fast to merge, others can be major pipeline updates.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="github_pr_response_time_plot" height="200"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> First response is when a comment is made by a GitHub user <em>other than</em> the original PR author</p>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="github_pr_response_time" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="github_pr_response_time" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>

    </div>
    <div class="col-lg-6">

        <h2 class="mt-0" id="github_issues">Issues<a href="#github_issues" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>GitHub issues can be created to log feature requests, bug reports or questions.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="github_issues_plot" height="200"></canvas>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="github_issues" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="github_issues" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>


        <h2 class="mt-0" id="github_issue_response_time">Issue response times<a href="#github_issue_response_time" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
        <p>A sign of an active community is a quick response time to issues. Here we see a frequency histogram of how long it takes to respond to and close issues.</p>
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="github_issue_response_time_plot" height="200"></canvas>
                <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> First response is when a comment is made by a GitHub user <em>other than</em> the original issue author</p>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="github_issue_response_time" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="github_issue_response_time" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>

    </div>
</div>

<h2 class="mt-0" id="contributor_leaderboard">Contributor Leaderboard<a href="#contributor_leaderboard" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
<p>We value each and every contribution to nf-core, no matter how small.
    However, that doesn't mean that we can't get competitive!</p>
<p>Here are the latest stats of who has contributed the greatest number of commits.
    The yellow bars show "core repositories" - repositories that are not pipelines
    (such as the code for this website!).
    A list of these repositories can be found <a href="#core_repos">below</a>.</p>
<div class="alert alert-light border small">
    <h6 class="alert-heading">Remember</h6>
    <ul>
        <li>There is more to contributing than commits! We're not counting issue comments, reviews or anything else here.</li>
        <li>People merging pull-requests get bonus commit counts from those merge commits.</li>
        <li>Some people commit often, others not so much. So it's not a perfect representation of amount of work - just a bit of fun!</li>
        <li><code>master</code> branch only, and all of the <a href="#caveats">other caveats</a>..</li>
    </ul>
</div>

<table class="table table-sm mt-5">
    <tbody>
        <?php
        $contributors = [];
        $contribution_counts = [];
        $contribution_counts_bytype = [];
        $top_repos = [];
        foreach (['pipelines', 'core_repos'] as $repo_type) {
            foreach ($stats_json->{$repo_type} as $repo_name => $repo) {
                foreach ($repo->contributors as $contributor) {
                    $login = $contributor->author->login;
                    $contributors[$login] = $contributor->author;
                    if (!isset($contribution_counts[$login])) {
                        $contribution_counts[$login] = 0;
                        $contribution_counts_bytype[$login] = ['pipelines' => 0, 'core_repos' => 0];
                        $top_repos[$login] = [$repo_name, $contributor->total];
                    }
                    $contribution_counts[$login] += $contributor->total;
                    $contribution_counts_bytype[$login][$repo_type] += $contributor->total;
                    if ($top_repos[$login][1] < $contributor->total) {
                        $top_repos[$login] = [$repo_name, $contributor->total];
                    }
                }
            }
        }
        arsort($contribution_counts);
        $max_count = max($contribution_counts);
        foreach ($contribution_counts as $login => $count) {
            $author = $contributors[$login];
            $pipeline_commits = $contribution_counts_bytype[$login]['pipelines'];
            $core_repo_commits = $contribution_counts_bytype[$login]['core_repos'];
            $pct_pipeline = ($pipeline_commits / $max_count) * 100;
            $pct_core = ($core_repo_commits / $max_count) * 100;
            echo '<tr>
      <td width="10%" class="pe-5">
        <a style="white-space: nowrap;" href="' .
                $author->html_url .
                '" target="_blank"><img src="' .
                $author->avatar_url .
                '" class="border rounded-circle me-1 mb-1" width="50" height="50"> @' .
                $author->login .
                '</a>
      </td>
      <td class="align-middle">
        <div class="progress" title="Pipelines: ' .
                $pipeline_commits .
                ' commits<br>Core repos: ' .
                $core_repo_commits .
                ' commits" data-bs-toggle="tooltip" data-html="true">
          <div class="progress-bar bg-success" role="progressbar" style="width: ' .
                $pct_pipeline .
                '%">' .
                ($pct_pipeline > 5 ? $pipeline_commits : '') .
                '</div>
          <div class="progress-bar bg-warning" role="progressbar" style="width: ' .
                $pct_core .
                '%">' .
                ($pct_core > 5 ? $core_repo_commits : '') .
                '</div>
        </div>
      </td>
      <td class="align-middle ps-5 small  font-monospace d-none d-md-table-cell" width="10%">
        <a href="/' .
                $top_repos[$login][0] .
                '" title="Repo with most commits (' .
                $top_repos[$login][1] .
                ' commits)" data-bs-toggle="tooltip">' .
                $top_repos[$login][0] .
                '
        <span class="badge rounded-pill bg-secondary float-end">' .
                $top_repos[$login][1] .
                '</span></a>
      </td>
    </tr>';
        }
        ?>
    </tbody>
</table>

<h2 class="mt-0" id="pipeline_numbers">Pipeline numbers<a href="#pipeline_numbers" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
<p>All nf-core pipelines are only considered stable when they have at least one release. Until then, they are classed as "in development".</p>

<div class="row">
    <div class="col-lg-6 offset-lg-3">
        <div class="card bg-body mt-4">
            <div class="card-body">
                <canvas id="pipeline_numbers_plot" height="200"></canvas>
                <p class="card-text small text-muted">
                    <a href="#" data-bs-target="pipeline_numbers" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-bs-target="pipeline_numbers" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
                </p>
            </div>
        </div>
    </div>
</div>

<?php // The pipeline and core repo tables are the same

foreach (['pipelines', 'core_repos'] as $repo_type): ?>

    <h2 class="mt-0" id="<?php echo $repo_type; ?>"><?php echo ucfirst(
    str_replace('_', ' ', $repo_type),
); ?><a href="#<?php echo $repo_type; ?>" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h2>
    <p class="text-info small">
        <i class="far fa-hand-point-right"></i>
        Click a row to see detailed statistics for that repository.
    </p>


    <div class="table-responsive">
        <table class="table table-hover table-sm small pipeline-stats-table">
            <thead class="">
                <tr>
                    <th>&nbsp;</th>
                    <th>Name</th>
                    <th>Age</th>
                    <?php if ($repo_type == 'pipelines'): ?><th class="">Releases</th><?php endif; ?>
                    <th class="">Committers</th>
                    <th class="">Commits</th>
                    <th class="">Stargazers</th>
                    <th class="">Watchers</th>
                    <th class="">Network Forks</th>
                    <th class="">Clones</th>
                    <th class="">Unique cloners</th>
                    <th class="">Repo views</th>
                    <th class="">Unique repo visitors</th>
                </tr>
            </thead>

            <tbody>
                <tr class="text-bold">
                    <th>&nbsp;</th>
                    <th>Total:</th>
                    <th class="font-weight-light"><?php echo count($pipelines); ?> pipelines</th>
                    <?php if ($repo_type == 'pipelines'): ?><th class="font-weight-light "><?php echo $stats_total[
    $repo_type
]['releases']; ?></th><?php endif; ?>
                    <th class="font-weight-light "><?php echo count(
                        $stats_total[$repo_type]['unique_committers'],
                    ); ?> unique</th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['total_commits']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['stargazers']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['watchers']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['network_forks_count']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['clones_count_total']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['clones_uniques_total']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['views_count_total']; ?></th>
                    <th class="font-weight-light "><?php echo $stats_total[$repo_type]['views_uniques_total']; ?></th>
                </tr>
                <?php echo implode($trows[$repo_type]); ?>
            </tbody>
            <tfoot class="">
                <tr>
                    <th>&nbsp;</th>
                    <th>Name</th>
                    <th>Age</th>
                    <?php if ($repo_type == 'pipelines'): ?><th class="">Releases</th><?php endif; ?>
                    <th class="">Committers</th>
                    <th class="">Commits</th>
                    <th class="">Stargazers</th>
                    <th class="">Watchers</th>
                    <th class="">Network Forks</th>
                    <th class="">Clones</th>
                    <th class="">Unique cloners</th>
                    <th class="">Repo views</th>
                    <th class="">Unique repo visitors</th>
                </tr>
            </tfoot>
        </table>
    </div>

<?php endforeach; ?>

<p class="mt-5 small text-muted">See also <a href="/pipeline_health">pipeline repository health</a>.</p>

<script type="text/javascript">
    var theme = "<?php echo $theme; ?>"
    $(function() {
        if (theme === 'auto') {
            if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
                theme = 'dark';
            } else {
                theme = 'light';
            }
        }
        // Placeholder for chart data
        var chartData = {};
        var charts = {};

        // Chart.JS base config
        var chartjs_base = {
            type: 'line',
            options: {
                title: {
                    display: true,
                    fontSize: 16
                },
                elements: {
                    line: {
                        borderWidth: 1,
                        tension: 0 // disables bezier curves
                    }
                },
                scales: {
                    xAxes: [{
                        type: 'time',
                        time: {
                            minUnit: 'day'
                        }
                    }],
                },
                legend: {
                    display: false
                },
                tooltips: {
                    mode: 'x'
                },
                plugins: {
                    zoom: {
                        zoom: {
                            enabled: true,
                            drag: true,
                            mode: 'x',
                            speed: 0.05
                        }
                    }
                }
            }
        };



        // Slack users chart
        chartData['slack'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['slack'].data = {
            datasets: [{
                    fill: 'origin', // explicitly fill the first dataset to the x axis
                    label: 'Inactive',
                    backgroundColor: 'rgba(150,150,150, 0.2)',
                    borderColor: 'rgba(150,150,150, 1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($stats_json->slack->user_counts as $timestamp => $users) {
                            echo '{ x: "' .
                                date('Y-m-d H:i:s', $timestamp) .
                                '", y: ' .
                                $users->inactive .
                                ' },' .
                                "\n\t\t\t";
                        } ?>
                    ]
                },
                {
                    label: 'Active',
                    backgroundColor: 'rgba(89, 37, 101, 0.2)',
                    borderColor: 'rgba(89, 37, 101, 1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($stats_json->slack->user_counts as $timestamp => $users) {
                            // Skip zeros (anything before 2010)
                            if ($timestamp < 1262304000) {
                                continue;
                            }
                            echo '{ x: "' .
                                date('Y-m-d H:i:s', $timestamp) .
                                '", y: ' .
                                $users->active .
                                ' },' .
                                "\n\t\t\t";
                        } ?>
                    ]
                }
            ]
        };
        chartData['slack'].options.title.text = 'nf-core Slack users over time';
        chartData['slack'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
        chartData['slack'].options.scales.yAxes = [{
            stacked: true
        }];
        chartData['slack'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('slack_users_plot').getContext('2d');
        charts['slack'] = new Chart(ctx, chartData['slack']);


        // GitHub org members chart
        chartData['gh_orgmembers'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['gh_orgmembers'].data = {
            datasets: [{
                backgroundColor: theme == 'light' ? 'rgba(0,0,0,0.2)' : 'rgba(255,255,255,0.7)',
                borderColor: theme == 'light' ? 'rgba(0,0,0,1)' : 'rgba(255,255,255,0.9)',
                pointRadius: 0,
                data: [
                    <?php foreach ($stats_json->gh_org_members as $timestamp => $count) {
                        // Skip zeros (anything before 2010)
                        if ($timestamp < 1262304000) {
                            continue;
                        }
                        echo '{ x: "' . date('Y-m-d H:i:s', $timestamp) . '", y: ' . $count . ' },' . "\n\t\t\t";
                    } ?>
                ]
            }]
        };
        chartData['gh_orgmembers'].options.title.text = 'nf-core GitHub organisation members over time';
        var ctx = document.getElementById('gh_orgmembers_plot').getContext('2d');
        charts['gh_orgmembers'] = new Chart(ctx, chartData['gh_orgmembers']);


        // GitHub contributors chart
        <?php
        asort($gh_contributors);
        $issues_cumulative_count = 0;
        $commits_cumulative_count = 0;
        $both_cumulative_count = 0;
        $contribs_commits = [];
        $contribs_issues = [];
        $contribs_both = [];
        foreach ($gh_contributors as $username => $timestamp) {
            // Make zeros and old timestamps start of 2018
            if ($timestamp < 1514764800) {
                $timestamp = 1514764800;
            }
            if (in_array($username, $gh_contributor_commits) && in_array($username, $gh_contributor_issues)) {
                $both_cumulative_count += 1;
            } elseif (in_array($username, $gh_contributor_commits)) {
                $commits_cumulative_count += 1;
            } elseif (in_array($username, $gh_contributor_issues)) {
                $issues_cumulative_count += 1;
            }
            $contribs_commits[] =
                '{ x: "' . date('Y-m-d H:i:s', $timestamp) . '", y: ' . $commits_cumulative_count . ' },' . "\n\t\t\t";
            $contribs_issues[] =
                '{ x: "' . date('Y-m-d H:i:s', $timestamp) . '", y: ' . $issues_cumulative_count . ' },' . "\n\t\t\t";
            $contribs_both[] =
                '{ x: "' . date('Y-m-d H:i:s', $timestamp) . '", y: ' . $both_cumulative_count . ' },' . "\n\t\t\t";
        }
        ?>
        chartData['gh_contribs'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['gh_contribs'].data = {
            datasets: [{
                    label: 'Commits',
                    backgroundColor: theme == 'light' ? 'rgba(0,0,0,0.2)' : 'rgba(255,255,255,0.7)',
                    borderColor: theme == 'light' ? 'rgba(0,0,0,1)' : 'rgba(255,255,255,0.9)',
                    pointRadius: 0,
                    fill: 'origin', // explicitly fill the first dataset to the x axis
                    data: [
                        <?php echo implode('', $contribs_commits); ?>
                    ]
                },
                {
                    label: 'Commits and Issues',
                    backgroundColor: 'rgba(104, 72, 186, 0.2)',
                    borderColor: 'rgba(104, 72, 186, 1.0)',
                    pointRadius: 0,
                    data: [
                        <?php echo implode('', $contribs_both); ?>
                    ]
                },
                {
                    label: 'Issues',
                    backgroundColor: 'rgba(83, 164, 81, 0.2)',
                    borderColor: 'rgba(83, 164, 81, 1.0)',
                    pointRadius: 0,
                    data: [
                        <?php echo implode('', $contribs_issues); ?>
                    ]
                }
            ]
        };
        chartData['gh_contribs'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
        chartData['gh_contribs'].options.scales.yAxes = [{
            stacked: true,
            scaleLabel: {
                display: true,
                labelString: 'Number of contributors'
            },
        }];
        chartData['gh_contribs'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        chartData['gh_contribs'].options.title.text = 'nf-core GitHub contributors over time';
        var ctx = document.getElementById('gh_contribs_plot').getContext('2d');
        charts['gh_contribs'] = new Chart(ctx, chartData['gh_contribs']);


        // Twitter followers chart
        chartData['twitter'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['twitter'].data = {
            datasets: [{
                backgroundColor: 'rgba(74, 161, 235, 0.2)',
                borderColor: 'rgba(74, 161, 235, 1)',
                pointRadius: 0,
                data: [
                    <?php foreach ($stats_json->twitter->followers_count as $timestamp => $count) {
                        // Skip zeros (anything before 2010)
                        if ($timestamp < 1262304000) {
                            continue;
                        }
                        echo '{ x: "' . date('Y-m-d H:i:s', $timestamp) . '", y: ' . $count . ' },' . "\n\t\t\t";
                    } ?>
                ]
            }]
        };
        chartData['twitter'].options.title.text = '@nf_core twitter followers users over time';
        var ctx = document.getElementById('twitter_followers_plot').getContext('2d');
        charts['twitter'] = new Chart(ctx, chartData['twitter']);

        // GitHub Pull Requests chart
        <?php
        $open_prs = [];
        $closed_prs = [];
        $dates_raw = array_unique(
            array_merge(
                array_keys($issues_json['stats']['prs']['daily_opened']),
                array_keys($issues_json['stats']['prs']['daily_closed']),
            ),
        );
        $dates = [];
        foreach ($dates_raw as $date) {
            $dates[strtotime($date)] = $date;
        }
        ksort($dates);
        $running_open_prs = 0;
        $running_closed_prs = 0;
        foreach ($dates as $ts => $date) {
            if (isset($issues_json['stats']['prs']['daily_opened'][$date])) {
                $running_open_prs += $issues_json['stats']['prs']['daily_opened'][$date];
            }
            if (isset($issues_json['stats']['prs']['daily_closed'][$date])) {
                $running_open_prs -= $issues_json['stats']['prs']['daily_closed'][$date];
                $running_closed_prs += $issues_json['stats']['prs']['daily_closed'][$date];
            }
            $open_prs[date('Y-m-d', $ts)] = $running_open_prs;
            $closed_prs[date('Y-m-d', $ts)] = $running_closed_prs;
        }
        ?>
        chartData['github_prs'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['github_prs'].data = {
            datasets: [{
                    label: 'Closed / Merged',
                    backgroundColor: 'rgba(104, 72, 186, 0.2)',
                    borderColor: 'rgba(104, 72, 186, 1)',
                    pointRadius: 0,
                    fill: 'origin', // explicitly fill the first dataset to the x axis
                    data: [
                        <?php foreach ($closed_prs as $date => $count) {
                            echo '{ x: "' . $date . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                },
                {
                    label: 'Open',
                    backgroundColor: 'rgba(83, 164, 81, 0.2)',
                    borderColor: 'rgba(83, 164, 81, 1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($open_prs as $date => $count) {
                            echo '{ x: "' . $date . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                }
            ]
        };
        chartData['github_prs'].options.title.text = 'GitHub Pull Requests over time';
        chartData['github_prs'].options.scales.yAxes = [{
            stacked: true
        }];
        chartData['github_prs'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
        chartData['github_prs'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('github_prs_plot').getContext('2d');
        charts['github_prs'] = new Chart(ctx, chartData['github_prs']);




        // GitHub issues chart
        <?php
        $open_issues = [];
        $closed_issues = [];
        // Get all dates for opening and closing issues
        $dates_raw = array_unique(
            array_merge(
                array_keys($issues_json['stats']['issues']['daily_opened']),
                array_keys($issues_json['stats']['issues']['daily_closed']),
            ),
        );
        // Date strings need sorting and formatting
        $dates = [];
        foreach ($dates_raw as $date) {
            $dates[strtotime($date)] = $date;
        }
        ksort($dates);
        $running_open_issues = 0;
        $running_closed_issues = 0;
        foreach ($dates as $ts => $date) {
            if (isset($issues_json['stats']['issues']['daily_opened'][$date])) {
                $running_open_issues += $issues_json['stats']['issues']['daily_opened'][$date];
            }
            if (isset($issues_json['stats']['issues']['daily_closed'][$date])) {
                $running_open_issues -= $issues_json['stats']['issues']['daily_closed'][$date];
                $running_closed_issues += $issues_json['stats']['issues']['daily_closed'][$date];
            }
            $open_issues[date('Y-m-d', $ts)] = $running_open_issues;
            $closed_issues[date('Y-m-d', $ts)] = $running_closed_issues;
        }
        ?>
        chartData['github_issues'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['github_issues'].data = {
            datasets: [{
                    label: 'Closed',
                    backgroundColor: 'rgba(199, 70, 78, 0.2)',
                    borderColor: 'rgba(199, 70, 78, 1)',
                    pointRadius: 0,
                    fill: 'origin', // explicitly fill the first dataset to the x axis
                    data: [
                        <?php foreach ($closed_issues as $date => $count) {
                            echo '{ x: "' . $date . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                },
                {
                    label: 'Open',
                    backgroundColor: 'rgba(83, 164, 81, 0.2)',
                    borderColor: 'rgba(83, 164, 81, 1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($open_issues as $date => $count) {
                            echo '{ x: "' . $date . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                }
            ]
        };
        chartData['github_issues'].options.title.text = 'GitHub Issues over time';
        chartData['github_issues'].options.scales.yAxes = [{
            stacked: true
        }];
        chartData['github_issues'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
        chartData['github_issues'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('github_issues_plot').getContext('2d');
        charts['github_issues'] = new Chart(ctx, chartData['github_issues']);




        <?php
        // Bin the data so that we can plot a histogram
        $hour = 60 * 60;
        $day = 60 * 60 * 24;
        $bins = [
            1 * $hour => '1 hour',
            2 * $hour => '2 hours',
            3 * $hour => '3 hours',
            4 * $hour => '4 hours',
            5 * $hour => '5 hours',
            6 * $hour => '6 hours',
            12 * $hour => '12 hours',
            1 * $day => '1 day',
            2 * $day => '2 days',
            3 * $day => '3 days',
            4 * $day => '4 days',
            5 * $day => '5 days',
            6 * $day => '6 days',
            7 * $day => '7 days',
            14 * $day => '14 days',
            31 * $day => '1 month',
            99999999999999999999 => 'Over 1 month',
        ];

        // ISSUES
        $gh_issue_response_hist = [];
        $gh_issue_close_hist = [];
        foreach ($issues_json['stats']['issues']['response_times'] as $rt) {
            // Find all bins that are bigger than the value, then take the smallest
            $key = min(
                array_filter(array_keys($bins), function ($ts) {
                    global $rt;
                    return $ts > $rt;
                }),
            );
            if (!isset($gh_issue_response_hist[$key])) {
                $gh_issue_response_hist[$key] = 0;
            }
            $gh_issue_response_hist[$key] += 1;
        }
        foreach ($issues_json['stats']['issues']['close_times'] as $rt) {
            // Find all bins that are bigger than the value, then take the smallest
            $key = min(
                array_filter(array_keys($bins), function ($ts) {
                    global $rt;
                    return $ts > $rt;
                }),
            );
            if (!isset($gh_issue_close_hist[$key])) {
                $gh_issue_close_hist[$key] = 0;
            }
            $gh_issue_close_hist[$key] += 1;
        }
        ksort($gh_issue_response_hist);
        ksort($gh_issue_close_hist);

        //PRs
        $gh_pr_response_hist = [];
        $gh_pr_close_hist = [];
        foreach ($issues_json['stats']['prs']['response_times'] as $rt) {
            // Find all bins that are bigger than the value, then take the smallest
            $key = min(
                array_filter(array_keys($bins), function ($ts) {
                    global $rt;
                    return $ts > $rt;
                }),
            );
            if (!isset($gh_pr_response_hist[$key])) {
                $gh_pr_response_hist[$key] = 0;
            }
            $gh_pr_response_hist[$key] += 1;
        }
        foreach ($issues_json['stats']['prs']['close_times'] as $rt) {
            // Find all bins that are bigger than the value, then take the smallest
            $key = min(
                array_filter(array_keys($bins), function ($ts) {
                    global $rt;
                    return $ts > $rt;
                }),
            );
            if (!isset($gh_pr_close_hist[$key])) {
                $gh_pr_close_hist[$key] = 0;
            }
            $gh_pr_close_hist[$key] += 1;
        }
        ksort($gh_pr_response_hist);
        ksort($gh_pr_close_hist);
        ?>
        // GitHub issues response time
        chartData['github_issue_response_time'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['github_issue_response_time'].data = {
            datasets: [{
                    label: 'Time to close',
                    backgroundColor: 'rgba(199, 70, 78, 0.5)',
                    borderColor: 'rgba(199, 70, 78, 0.1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($bins as $key => $label) {
                            echo ($gh_issue_close_hist[$key] / count($issues_json['stats']['issues']['close_times'])) *
                                100 .
                                ', ';
                        } ?>
                    ]
                },
                {
                    label: 'Time to first response',
                    backgroundColor: 'rgba(83, 164, 81, 0.5)',
                    borderColor: 'rgba(83, 164, 81, 0.1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($bins as $key => $label) {
                            echo ($gh_issue_response_hist[$key] /
                                count($issues_json['stats']['issues']['response_times'])) *
                                100 .
                                ', ';
                        } ?>
                    ]
                }
            ]
        };
        chartData['github_issue_response_time'].type = 'bar';
        chartData['github_issue_response_time'].options.scales.xAxes = [{
            type: 'category',
            labels: ["<?php echo implode('", "', array_values($bins)); ?>"],
        }];
        chartData['github_issue_response_time'].options.scales.yAxes = [{
            scaleLabel: {
                display: true,
                labelString: 'Percentage of issues'
            },
            ticks: {
                // Include a dollar sign in the ticks
                callback: function(value, index, values) {
                    return value + '%';
                }
            }
        }];
        // Round the tooltip values
        chartData['github_issue_response_time'].options.tooltips = {
            mode: 'index',
            callbacks: {
                label: function(tooltipItem, data) {
                    var label = data.datasets[tooltipItem.datasetIndex].label || '';
                    if (label) {
                        label += ': ';
                    }
                    label += Math.round(tooltipItem.yLabel * 100) / 100;
                    label += '%';
                    return label;
                }
            }
        }
        chartData['github_issue_response_time'].options.plugins.zoom = false;
        chartData['github_issue_response_time'].options.title.text = 'GitHub Issues Response Time';
        chartData['github_issue_response_time'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('github_issue_response_time_plot').getContext('2d');
        charts['github_issue_response_time'] = new Chart(ctx, chartData['github_issue_response_time']);




        // GitHub PRs response time
        chartData['github_pr_response_time'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['github_pr_response_time'].data = {
            datasets: [{
                    label: 'Time to merge / close',
                    backgroundColor: 'rgba(104, 72, 186, 0.5)',
                    borderColor: 'rgba(104, 72, 186, 0.1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($bins as $key => $label) {
                            echo ($gh_pr_close_hist[$key] / count($issues_json['stats']['prs']['close_times'])) * 100 .
                                ', ';
                        } ?>
                    ]
                },
                {
                    label: 'Time to first response',
                    backgroundColor: 'rgba(83, 164, 81, 0.5)',
                    borderColor: 'rgba(83, 164, 81, 0.1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($bins as $key => $label) {
                            echo ($gh_pr_response_hist[$key] / count($issues_json['stats']['prs']['response_times'])) *
                                100 .
                                ', ';
                        } ?>
                    ]
                }
            ]
        };
        chartData['github_pr_response_time'].type = 'bar';
        chartData['github_pr_response_time'].options.scales.xAxes = [{
            type: 'category',
            labels: ["<?php echo implode('", "', array_values($bins)); ?>"],
        }];
        chartData['github_pr_response_time'].options.scales.yAxes = [{
            scaleLabel: {
                display: true,
                labelString: 'Percentage of PRs'
            },
            ticks: {
                // Include a dollar sign in the ticks
                callback: function(value, index, values) {
                    return value + '%';
                }
            }
        }];
        // Round the tooltip values
        chartData['github_pr_response_time'].options.tooltips = {
            mode: 'index',
            callbacks: {
                label: function(tooltipItem, data) {
                    var label = data.datasets[tooltipItem.datasetIndex].label || '';
                    if (label) {
                        label += ': ';
                    }
                    label += Math.round(tooltipItem.yLabel * 100) / 100;
                    label += '%';
                    return label;
                }
            }
        }
        chartData['github_pr_response_time'].options.plugins.zoom = false;
        chartData['github_pr_response_time'].options.title.text = 'GitHub Pull Request Response Time';
        chartData['github_pr_response_time'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('github_pr_response_time_plot').getContext('2d');
        charts['github_pr_response_time'] = new Chart(ctx, chartData['github_pr_response_time']);


        // Pipeline numbers plot
        <?php
        $pipeline_dates = [];
        foreach ($pipelines as $pipeline) {
            $pipeline_dates[strtotime($pipeline->created_at)] = 'create';
            $first_release = false;
            foreach ($pipeline->releases as $release) {
                if (!$first_release) {
                    $first_release = strtotime($release->published_at);
                }
                $first_release = min($first_release, strtotime($release->published_at));
            }
            if ($first_release) {
                $pipeline_dates[$first_release] = 'release';
            }
        }
        ksort($pipeline_dates);
        $dev_pipelines = [];
        $released_pipelines = [];
        $dev_running = 0;
        $release_running = 0;
        foreach ($pipeline_dates as $date => $dtype) {
            if ($dtype == 'create') {
                $dev_running += 1;
            }
            if ($dtype == 'release') {
                $dev_running -= 1;
                $release_running += 1;
            }
            $released_pipelines[$date] = $release_running;
            $dev_pipelines[$date] = $dev_running;
        }
        ?>
        chartData['pipeline_numbers'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['pipeline_numbers'].data = {
            datasets: [{
                    label: 'Released',
                    backgroundColor: 'rgba(83, 164, 81, 0.3)',
                    borderColor: 'rgba(83, 164, 81, 1)',
                    pointRadius: 0,
                    fill: 'origin', // explicitly fill the first dataset to the x axis
                    data: [
                        <?php foreach ($released_pipelines as $timestamp => $count) {
                            echo '{ x: "' . date('Y-m-d', $timestamp) . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                },
                {
                    label: 'In development',
                    backgroundColor: 'rgba(199, 70, 78, 0.1)',
                    borderColor: 'rgba(199, 70, 78, 1)',
                    pointRadius: 0,
                    data: [
                        <?php foreach ($dev_pipelines as $timestamp => $count) {
                            echo '{ x: "' . date('Y-m-d', $timestamp) . '", y: ' . $count . ' },' . "\n\t\t\t";
                        } ?>
                    ]
                }
            ]
        };
        chartData['pipeline_numbers'].options.title.text = 'nf-core pipeline numbers over time';
        chartData['pipeline_numbers'].options.scales.yAxes = [{
            stacked: true
        }];
        chartData['pipeline_numbers'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
        chartData['pipeline_numbers'].options.legend = {
            position: 'bottom',
            labels: {
                lineWidth: 1
            }
        };
        var ctx = document.getElementById('pipeline_numbers_plot').getContext('2d');
        charts['pipeline_numbers'] = new Chart(ctx, chartData['pipeline_numbers']);



        // Repo clones plot
        chartData['repo_clones'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['repo_clones'].data = {
            datasets: [{
                    label: 'Clones per day',
                    backgroundColor: 'rgba(83, 164, 81, 1.0)',
                    borderColor: 'rgba(83, 164, 81, 1.0)',
                    fill: false,
                    pointRadius: 0,
                    yAxisID: 'y-axis-count',
                    data: [
                        <?php
                        foreach( $clone_counts as $idx => $count ){
                            echo '{ x: "' . date('Y-m-d', strtotime($count['timestamp'])) . '", y: ' . $count['sum_total_clones']  . ' },' . "\n\t\t\t";
                        }
                        ?>
                    ]
                },
                {
                    label: 'Unique cloners per day',
                    backgroundColor: 'rgba(33, 94, 190, 1.0)',
                    borderColor: 'rgba(33, 94, 190, 1.0)',
                    fill: false,
                    pointRadius: 0,
                    yAxisID: 'y-axis-uniques',
                    data: [
                        <?php
                        foreach( $clone_counts as $idx => $count ){
                            echo '{ x: "' . date('Y-m-d', strtotime($count['timestamp'])) . '", y: ' . $count['sum_total_clones_unique']  . ' },' . "\n\t\t\t";
                        }
                        ?>
                    ]
                }
            ]
        };
        chartData['repo_clones'].options.title.text = 'nf-core git clones per day';
        chartData['repo_clones'].options.elements.line.borderWidth = 2;
        chartData['repo_clones'].options.scales.yAxes = [{
                id: 'y-axis-count',
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Clones per day',
                    fontColor: 'rgba(83, 164, 81, 1.0)',
                },
                position: 'left',
                gridLines: {
                    drawOnChartArea: false
                }
            },
            {
                id: 'y-axis-uniques',
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Unique cloners per day',
                    fontColor: 'rgba(33, 94, 190, 1.0)',
                },
                position: 'right',
                gridLines: {
                    drawOnChartArea: false
                }
            }
        ];
        var ctx = document.getElementById('repo_clones_plot').getContext('2d');
        charts['repo_clones'] = new Chart(ctx, chartData['repo_clones']);



        // Repo views plot
        chartData['repo_views'] = JSON.parse(JSON.stringify(chartjs_base));
        chartData['repo_views'].data = {
            datasets: [{
                    label: 'Views per day',
                    backgroundColor: 'rgba(83, 164, 81, 1.0)',
                    borderColor: 'rgba(83, 164, 81, 1.0)',
                    fill: false,
                    pointRadius: 0,
                    yAxisID: 'y-axis-count',
                    data: [
                        <?php
                        // add data point for each timestamp in view_count
                        foreach( $view_counts as $idx => $count ){
                            echo '{ x: "' . date('Y-m-d', strtotime($count['timestamp'])) . '", y: ' . $count['sum_total_views']  . ' },' . "\n\t\t\t";
                        }
                        // ksort($stats_total_allrepos['views_count']);
                        // foreach ($stats_total_allrepos['views_count'] as $timestamp => $count) {
                        //     $timestamp = strtotime($timestamp);
                        //     // Skip zeros (anything before 2010)
                        //     if ($timestamp < 1262304000) {
                        //         continue;
                        //     }
                        //     echo '{ x: "' . date('Y-m-d', $timestamp) . '", y: ' . $count . ' },' . "\n\t\t\t";
                        // }
                        ?>
                    ]
                },
                {
                    label: 'Unique visitors per day',
                    backgroundColor: 'rgba(33, 94, 190, 1.0)',
                    borderColor: 'rgba(33, 94, 190, 1.0)',
                    fill: false,
                    pointRadius: 0,
                    yAxisID: 'y-axis-uniques',
                    data: [
                        <?php
                        foreach( $view_counts as $idx => $count ){
                            echo '{ x: "' . date('Y-m-d', strtotime($count['timestamp'])) . '", y: ' . $count['sum_total_views_unique']  . ' },' . "\n\t\t\t";
                        }
                        ?>
                    ]
                }
            ]
        };
        chartData['repo_views'].options.title.text = 'nf-core repository web views per day';
        chartData['repo_views'].options.elements.line.borderWidth = 2;
        chartData['repo_views'].options.scales.yAxes = [{
                id: 'y-axis-count',
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Views per day',
                    fontColor: 'rgba(83, 164, 81, 1.0)',
                },
                position: 'left',
                gridLines: {
                    drawOnChartArea: false
                }
            },
            {
                id: 'y-axis-uniques',
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Unique visitors per day',
                    fontColor: 'rgba(33, 94, 190, 1.0)',
                },
                position: 'right',
                gridLines: {
                    drawOnChartArea: false
                }
            }
        ];
        var ctx = document.getElementById('repo_views_plot').getContext('2d');
        charts['repo_views'] = new Chart(ctx, chartData['repo_views']);

        // Make canvas2svg work with ChartJS
        // https://stackoverflow.com/a/52151467/713980
        function canvas2svgTweakLib() {
            C2S.prototype.getContext = function(contextId) {
                if (contextId == "2d" || contextId == "2D") {
                    return this;
                }
                return null;
            }
            C2S.prototype.style = function() {
                return this.__canvas.style
            }
            C2S.prototype.getAttribute = function(name) {
                return this[name];
            }
            C2S.prototype.addEventListener = function(type, listener, eventListenerOptions) {
                console.log("canvas2svg.addEventListener() not implemented.")
            }
        }
        canvas2svgTweakLib();

        function exportChartJsSVG(target) {

            // Bump default font size
            Chart.defaults.global.defaultFontSize = 18;

            // Turn off responsiveness
            chartData[target].options.responsive = false;
            chartData[target].options.animation = false;
            chartData[target].options.plugins.zoom = false;
            chartData[target].options.plugins.zoom = false;

            // Add extra height if we have a legend
            canvas_height = 400;
            if (chartData[target].options.legend.position == 'bottom') {
                canvas_height = 450;
            }
            // Made plot wider if it's views or clones
            canvas_width = 800;
            if (chartData[target].options.elements.line.borderWidth == 2) {
                canvas_width = 1500;
            }

            // canvas2svg 'mock' context
            var svgContext = C2S(canvas_width, canvas_height);
            // new chart on 'mock' context fails:
            var mySvg = new Chart(svgContext, chartData[target]);
            // Failed to create chart: can't acquire context from the given item
            var svg = svgContext.getSerializedSvg(true);
            // Trigger browser download with SVG
            var blob = new Blob([svg], {
                type: "text/plain;charset=utf-8"
            });
            saveAs(blob, 'nf-core_' + target + '_plot.svg');

            // Reset plots to defaults
            Chart.defaults.global.defaultFontSize = 12;
            chartData[target].options.responsive = true;
            chartData[target].options.animation = true;
            chartData[target].options.plugins.zoom = {
                zoom: {
                    enabled: true,
                    drag: true,
                    mode: 'x',
                    speed: 0.05
                }
            }
        }
        $('.dl_plot_svg').click(function(e) {
            e.preventDefault();
            var target = $(this).data('bsTarget');
            exportChartJsSVG(target);
        });
        $('.reset_chart_zoom').click(function(e) {
            e.preventDefault();
            var target = $(this).data('bsTarget');
            charts[target].resetZoom();
        });

    });
</script>

<?php
// Stop the clocks!
$end_time = microtime(true);
$time_taken = round($end_time - $start_time, 5);

$subfooter =
    '<p class="mb-0"><i class="far fa-clock"></i> Last updated: ' .
    date('d-m-Y', $stats_json->updated) .
    '. Page generated in ' .
    $time_taken .
    ' seconds.</p>';

// include '../includes/footer.php';

file_put_contents('stats_static.html', ob_get_contents());


?>
