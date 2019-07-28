<?php
// ini_set('display_errors', 1);
// ini_set('display_startup_errors', 1);
// error_reporting(E_ALL);

function round_nicely($num){
  if($num > 1000000){
    $num /= 1000000;
    $num = round($num, 2).'M';
  } else if($num > 1000){
    $num /= 1000;
    $num = round($num, 2).'K';
  }
  return $num;
}

$title = 'nf-core in numbers';
$subtitle = 'Measuring activity across the nf-core community.';
$import_chartjs = true;
include('../includes/header.php');

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

$stats_json_fn = dirname(dirname(__FILE__)).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn));

// Convenience variables
$slack_users = $stats_json->slack->user_counts->{$stats_json->updated};
$twitter_datekeys = array_keys(get_object_vars($stats_json->twitter->followers_count));
$twitter_users = $stats_json->twitter->followers_count->{max($twitter_datekeys)};

# echo '<pre>'.print_r($stats, true).'</pre>';

// This is horrible, I know sorry. I'm in a rush :(
$stats_total_allrepos = [
  'clones_count' => array(),
  'clones_uniques' => array(),
  'views_count' => array(),
  'views_uniques' => array()
];

// Run everything twice - keep pipelines and core repos seperate
foreach(['pipelines', 'core_repos'] as $repo_type):

$stats = $stats_json->{$repo_type};
$stats_total[$repo_type] = [
  'releases' => 0,
  'stargazers' => 0,
  'watchers' => 0,
  'forks' => 0,
  'clones_count' => array(),
  'clones_uniques' => array(),
  'views_count' => array(),
  'views_uniques' => array(),
  'clones_count_total' => 0,
  'clones_uniques_total' => 0,
  'views_count_total' => 0,
  'views_uniques_total' => 0,
  'unique_contributors' => [],
  'total_commits' => 0,
];

$trows[$repo_type] = [];
foreach($stats as $repo_name => $repo):
  // Exit quietly if something has gone wrong
  if(!isset($repo->repo_metrics)){
    echo '<!-- ERROR: $repo->repo_metrics not set for "'.$repo_name.'" -->';
    continue;
  }
  if(!isset($repo->repo_metrics->{$stats_json->updated})){
    echo '<!-- ERROR: $repo->repo_metrics->'.$stats_json->updated.' not set for "'.$repo_name.'" -->';
    continue;
  }
  // Ok, continue!
  $metrics = $repo->repo_metrics->{$stats_json->updated};
  $stats_total[$repo_type]['releases'] += isset($repo->num_releases) ? $repo->num_releases : 0;
  $stats_total[$repo_type]['stargazers'] += $metrics->stargazers_count;
  $stats_total[$repo_type]['forks'] += $metrics->forks_count;

  foreach(['clones_count', 'clones_uniques', 'views_count', 'views_uniques'] as $key){
    if(isset($repo->{$key})){
      foreach($repo->{$key} as $timestamp => $count){
        if(!isset($stats_total_allrepos[$key][$timestamp])){
          $stats_total_allrepos[$key][$timestamp] = 0;
        }
        $stats_total_allrepos[$key][$timestamp] += $count;
        if(!isset($stats_total[$repo_type][$key][$timestamp])){
          $stats_total[$repo_type][$key][$timestamp] = 0;
        }
        $stats_total[$repo_type][$key][$timestamp] += $count;
      }
    }
    $stats_total[$repo_type][$key.'_total'] += $repo->{$key.'_total'};
  }

  $total_commits = 0;
  foreach($repo->contributors as $contributor){
    $gh_username = $contributor->author->login;
    $stats_total[$repo_type]['unique_contributors'][$gh_username] = 0;
    $stats_total['total']['unique_contributors'][$gh_username] = 0;
    $stats_total[$repo_type]['total_commits'] += $contributor->total;
    $total_commits += $contributor->total;
  }
  ob_start();
  ?>
  <tr>
    <td><?php
    if($metrics->archived){
      echo '<small class="status-icon text-warning ml-2 fas fa-archive" data-toggle="tooltip" aria-hidden="true" title="This repo has been archived and is no longer being maintained."></small>';
    } else if($repo_type == 'pipelines'){
      // Edge case where a new pipeline is added but stats hasn't rerun yet
      if(!isset($repo->num_releases)){ $repo->num_releases = 0; }
      if($repo->num_releases){
        echo '<small class="status-icon text-success ml-2 fas fa-check" data-toggle="tooltip" aria-hidden="true" title="This pipeline is released, tested and good to go."></small>';
      } else {
        echo '<small class="status-icon text-danger ml-2 fas fa-wrench" data-toggle="tooltip" aria-hidden="true" title="This pipeline is under active development. Once released on GitHub, it will be production-ready."></small>';
      }
    }
    ?></td>
    <td><?php echo '<a href="'.$metrics->html_url.'" target="_blank"><span class="d-none d-lg-inline">nf-core/</span>'.$metrics->name.'</a>'; ?></td>
    <td><?php echo time_ago($metrics->created_at, false); ?></td>
    <?php if($repo_type == 'pipelines'): ?><td class="text-right"><?php echo $repo->num_releases; ?></td><?php endif; ?>
    <td class="text-right"><?php echo $repo->num_contributors; ?></td>
    <td class="text-right"><?php echo $total_commits; ?></td>
    <td class="text-right"><?php echo $metrics->stargazers_count; ?></td>
    <td class="text-right"><?php echo $metrics->forks_count; ?></td>
    <td class="text-right"><?php echo $repo->clones_count_total; ?></td>
    <td class="text-right"><?php echo $repo->clones_uniques_total; ?></td>
    <td class="text-right"><?php echo $repo->views_count_total; ?></td>
    <td class="text-right"><?php echo $repo->views_uniques_total; ?></td>
  </tr>
<?php
$trows[$repo_type][] = ob_get_contents();
ob_end_clean();
endforeach;

endforeach;

foreach(array_keys($stats_total['pipelines']) as $akey){
  if($akey == 'unique_contributors'){
    continue;
  }
  $stats_total['total'][$akey] = $stats_total['pipelines'][$akey] + $stats_total['core_repos'][$akey];
}

//
//
// TOP CARD DECK
//
//
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
      <li><a href="#pipelines">Pipelines</a></li>
      <li><a href="#core_repos">Core repositories</a></li>
    </ul>
  </li>
</ul>

<section id="caveats">
  <div class="card mb-3 mt-2">
    <div class="card-header">
      <a href="#caveats_<?php echo $repo_type; ?>" data-toggle="collapse" data-target="#caveats_<?php echo $repo_type; ?>" class="text-muted">
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
          <li>Metrics are fetched once per day (last checked <?php echo date('Y-m-d', $stats_json->updated); ?>).</li>
        </ul>
      </div>
    </div>
  </div>
</section>

<section id="community">
<h1>Community</h1>
<p>The numbers below track our growth over the various channels that the nf-core community operates in.</p>
<p class="text-info small">
  <i class="far fa-hand-point-right"></i>
  Click a number to see how the community has grown over time
</p>

<div class="card-group text-center stats_keynumbers">
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#slack" class="text-body text-decoration-none stretched-link"><?php echo $slack_users->total; ?></a></p>
      <p class="card-text text-muted">Slack users</p>
    </div>
    <div class="bg-icon" style="color: rgba(89, 37, 101, 0.1);"><i class="fab fa-slack"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#gh_orgmembers" class="text-body text-decoration-none stretched-link"><?php echo $stats_json->gh_org_members->{$stats_json->updated}; ?></a></p>
      <p class="card-text text-muted">GitHub organisation members</p>
    </div>
    <div class="bg-icon"><i class="fab fa-github"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#gh_contribs" class="text-body text-decoration-none stretched-link"><?php echo count($stats_total['total']['unique_contributors']); ?></a></p>
      <p class="card-text text-muted">GitHub contributors</p>
    </div>
    <div class="bg-icon"><i class="fas fa-code-branch"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#twitter" class="text-body text-decoration-none stretched-link"><?php echo $twitter_users; ?></a></p>
      <p class="card-text text-muted">Twitter followers</p>
    </div>
    <div class="bg-icon" style="color: rgba(74, 161, 235, 0.2);"><i class="fab fa-twitter"></i></div>
  </div>
</div>

<div class="row">
  <div class="col-lg-6">

    <section id="slack">
      <h2>Slack</h2>
      <p>Slack is a real-time messaging tool, with discussion split into channels and groups.
      We use it to provide help to people running nf-core pipelines, as well as discussing development ideas.
      You can join the nf-core slack by getting an invite <a href="https://nf-core-invite.herokuapp.com/">here</a>.</p>
      <div class="card bg-light mt-4">
        <div class="card-body">
          <canvas id="slack_users_plot" height="200"></canvas>
          <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> Slack considers users to be inactive when they haven't used slack for the previous 14 days.</p>
          <p class="card-text small text-muted mb-1"><i class="fas fa-exclamation-triangle"></i> Data from before 2019-07-24 fudged by reverse-engineering billing details on the slack admin pages.</p>
          <p class="card-text small text-muted">
            <a href="#" data-target="slack" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-target="slack" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
          </p>
        </div>
      </div>
    </section> <!-- <section id="slack"> -->

    <section id="twitter">
      <h2>Twitter followers</h2>
      <p>We use our <a href="https://twitter.com/nf_core">@nf_core</a> twitter account to send automated tweets about new pipeline releases and other updates relevant to the community.
      Follower counts give some indication to the level of interest in the nf-core project.</p>
      <div class="card bg-light mt-4">
        <div class="card-body">
          <canvas id="twitter_followers_plot" height="150"></canvas>
          <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-exclamation-triangle"></i> Data from before 2019-06-26 fudged by reverse-engineering a tiny sparkline plot on the twitter analytics website.</p>
          <p class="card-text small text-muted"><a href="#" data-target="twitter" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a>  &nbsp;/&nbsp; <a href="#" data-target="twitter" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
        </div>
      </div>
    </section> <!-- <section id="twitter"> -->
  </div>

  <div class="col-lg-6">
    <section id="gh_orgmembers">

      <h2>GitHub organisation members</h2>
      <p>We use GitHub to manage all of the code written for nf-core.
      It's a fantastic platform and provides a huge number of tools.
      We have a GitHub organisation called <a href="https://github.com/nf-core/">nf-core</a> which anyone can join:
      drop us a note <a href="https://github.com/nf-core/nf-co.re/issues/3">here</a> or anywhere and we'll send you an invite.
      </p>
      <p>It's not required to be a member of the nf-core GitHub organisation to contribute.
      However, members get the nf-core logo listed on their profile page and full write-access to all nf-core repositories.
      </p>
      <div class="card bg-light mt-4">
        <div class="card-body">
          <canvas id="gh_orgmembers_plot" height="150"></canvas>
          <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-exclamation-triangle"></i> By default, organisation membership is private. This is why you'll see a lower number if you visit the <a href="https://github.com/nf-core/">nf-core organisation page</a> and are not a member.
          <p class="card-text small text-muted"><a href="#" data-target="gh_orgmembers" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-target="gh_orgmembers" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
        </div>
      </div>
    </section> <!-- <section id="gh_orgmembers"> -->

    <section id="gh_contribs">
      <h2>GitHub contributors</h2>
      <p>Anybody can fork nf-core repositories and open a pull-request.
      Here we count how many different people have contributed at least one commit to an nf-core repository.</p>
      <div class="card bg-light mt-4">
        <div class="card-body">
          <canvas id="gh_contribs_plot" height="150"></canvas>
          <p class="card-text small text-muted mt-3 mb-1"><i class="fas fa-info-circle"></i> Some pipelines have been moved to the nf-core organisation instead of being forked. Contributions for these repos may predate nf-core.</p>
          <p class="card-text small text-muted"><a href="#" data-target="gh_contribs" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> Download as SVG</a> &nbsp;/&nbsp; <a href="#" data-target="gh_contribs" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a></p>
        </div>
      </div>
    </section> <!-- <section id="gh_contribs"> -->

  </div>
</div>

</section> <!-- <section id="community"> -->


<section id="code">
<h1>Code stats</h1>
<p>Whilst we always prefer quality over quantity, these numbers reflect the work output from the nf-core community.</p>

<div class="card-group text-center stats_keynumbers">
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><?php echo count(get_object_vars($stats_json->pipelines)) + count(get_object_vars($stats_json->core_repos)); ?></p>
      <p class="card-text text-muted">Repositories</p>
    </div>
    <div class="bg-icon"><i class="far fa-folder"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $stats_total['pipelines']['releases']; ?></p>
      <p class="card-text text-muted">Pipeline releases</p>
    </div>
    <div class="bg-icon"><i class="fas fa-tags"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><?php echo round_nicely($stats_total['pipelines']['total_commits'] + $stats_total['core_repos']['total_commits']); ?></p>
      <p class="card-text text-muted">Commits</p>
    </div>
    <div class="bg-icon"><i class="far fa-file-code"></i></div>
  </div>
</div>

<section id="repo_traffic">
  <h2>Repository traffic</h2>
  <p>Every time a nextflow user pulls an nf-core pipeline, the repository is cloned. Here we can track how much that happens across all nf-core repositories.
  Please note that these numbers come with some caveats <a href="#caveats">[ see more ]</a>.</p>
  <p>Additionally, GitHub tracks how many times people view repository web pages on github.com.</p>

  <div class="card mt-4">
    <div class="card-header">
      <span class="float-right small text-muted">
        <a href="#" data-target="repo_clones" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
        &nbsp;/&nbsp; <a href="#" data-target="repo_clones" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
      </span>
      Git clones: All nf-core repositories
    </div>
    <div class="card-body">
      <canvas id="repo_clones_plot" height="80"></canvas>
    </div>
    <div class="card-footer text-muted text-center small">
      <div class="row">
        <div class="col-6 border-right border-secondary">
          <span class="text-body lead"><?php echo round_nicely($stats_total['pipelines']['clones_count_total'] + $stats_total['core_repos']['clones_count_total']); ?></span><br>Clones
        </div>
        <div class="col-6" data-toggle="tooltip" title="Note: Unique per repository. Will double-count the same person cloning two different repositories.">
          <span class="text-body lead"><?php echo round_nicely($stats_total['pipelines']['clones_uniques_total'] + $stats_total['core_repos']['clones_uniques_total']); ?></span><br>Unique cloners
        </div>
      </div>
    </div>
  </div>

  <div class="card mt-4">
    <div class="card-header">
      <span class="float-right small text-muted">
        <a href="#" data-target="repo_views" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
        &nbsp;/&nbsp; <a href="#" data-target="repo_views" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
      </span>
      Visitors: All nf-core repositories
    </div>
    <div class="card-body">
      <canvas id="repo_views_plot" height="80"></canvas>
    </div>
    <div class="card-footer text-muted text-center small">
      <div class="row">
        <div class="col-6 border-right border-secondary">
          <span class="text-body lead"><?php echo round_nicely($stats_total['pipelines']['views_count_total'] + $stats_total['core_repos']['views_count_total']); ?></span><br>Views
        </div>
        <div class="col-6" data-toggle="tooltip" title="Note: Unique per repository. Will double-count the same person viewing two different repositories.">
          <span class="text-body lead"><?php echo round_nicely($stats_total['pipelines']['views_uniques_total'] + $stats_total['core_repos']['views_uniques_total']); ?></span><br>Unique visitors
        </div>
      </div>
    </div>
  </div>
</section> <!-- <section id="repo_views"> -->

<?php
// The pipeline and core repo tables are the same
foreach(['pipelines', 'core_repos'] as $repo_type): ?>

<section id="<?php echo $repo_type; ?>">
<h2><?php echo ucfirst(str_replace('_', ' ', $repo_type)); ?></h2>
<p class="text-info small">
  <i class="far fa-hand-point-right"></i>
  Click a row to see detailed statistics for that repository.
</p>


<div class="table-responsive">
  <table class="table table-hover table-sm small pipeline-stats-table">
    <thead class="thead-light">
      <tr>
        <th>&nbsp;</th>
        <th>Name</th>
        <th>Age</th>
        <?php if($repo_type == 'pipelines'): ?><th class="text-right">Releases</th><?php endif; ?>
        <th class="text-right">Contributors</th>
        <th class="text-right">Commits</th>
        <th class="text-right">Stargazers</th>
        <th class="text-right">Forks</th>
        <th class="text-right">Clones</th>
        <th class="text-right">Unique cloners</th>
        <th class="text-right">Repo views</th>
        <th class="text-right">Unique repo visitors</th>
      </tr>
    </thead>
    <thead class="thead-dark">
      <tr>
        <th>&nbsp;</th>
        <th>Total:</th>
        <th class="font-weight-light"><?php echo count($pipelines); ?> pipelines</th>
        <?php if($repo_type == 'pipelines'): ?><th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['releases']; ?></th><?php endif; ?>
        <th class="font-weight-light text-right"><?php echo count($stats_total[$repo_type]['unique_contributors']); ?> unique</th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['total_commits']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['stargazers']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['forks']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['clones_count_total']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['clones_uniques_total']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['views_count_total']; ?></th>
        <th class="font-weight-light text-right"><?php echo $stats_total[$repo_type]['views_uniques_total']; ?></th>
      </tr>
    </thead>
    <tbody>
    <?php echo implode($trows[$repo_type]); ?>
    </tbody>
    <tfoot class="thead-light">
      <tr>
        <th>&nbsp;</th>
        <th>Name</th>
        <th>Age</th>
        <?php if($repo_type == 'pipelines'): ?><th class="text-right">Releases</th><?php endif; ?>
        <th class="text-right">Contributors</th>
        <th class="text-right">Commits</th>
        <th class="text-right">Stargazers</th>
        <th class="text-right">Forks</th>
        <th class="text-right">Clones</th>
        <th class="text-right">Unique cloners</th>
        <th class="text-right">Repo views</th>
        <th class="text-right">Unique repo visitors</th>
      </tr>
    </tfoot>
  </table>
</div>
</section> <!-- <section id="<?php echo $repo_type; ?>"> -->

<?php endforeach; ?>

</section> <!-- <section id="code"> -->

<script type="text/javascript">
$(function(){

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
        xAxes: [{ type: 'time' }]
      },
      legend: {
        display: false
      },
      tooltips: { mode: 'x' },
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
    datasets: [
      {
        label: 'Inactive',
        backgroundColor: 'rgba(150,150,150, 0.2)',
        borderColor: 'rgba(150,150,150, 1)',
        pointRadius: 0,
        data: [
          <?php
          foreach($stats_json->slack->user_counts as $timestamp => $users){
            echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$users->inactive.' },'."\n\t\t\t";
          }
          ?>
        ]
      },
      {
        label: 'Active',
        backgroundColor: 'rgba(89, 37, 101, 0.2)',
        borderColor: 'rgba(89, 37, 101, 1)',
        pointRadius: 0,
        data: [
          <?php
          foreach($stats_json->slack->user_counts as $timestamp => $users){
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$users->active.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['slack'].options.title.text = 'nf-core Slack users over time';
  chartData['slack'].options.scales.yAxes = [{stacked: true }];
  chartData['slack'].options.legend = {
    position: 'bottom',
    labels: { lineWidth: 1 }
  };
  var ctx = document.getElementById('slack_users_plot').getContext('2d');
  charts['slack'] = new Chart(ctx, chartData['slack']);


  // GitHub org members chart
  chartData['gh_orgmembers'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['gh_orgmembers'].data = {
    datasets: [
      {
        backgroundColor: 'rgba(0,0,0,0.2)',
        borderColor: 'rgba(0,0,0,1)',
        pointRadius: 0,
        data: [
          <?php
          foreach($stats_json->gh_org_members as $timestamp => $count){
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['gh_orgmembers'].options.title.text = 'nf-core GitHub organisation members over time';
  var ctx = document.getElementById('gh_orgmembers_plot').getContext('2d');
  charts['gh_orgmembers'] = new Chart(ctx, chartData['gh_orgmembers']);


  // GitHub contributors chart
  chartData['gh_contribs'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['gh_contribs'].data = {
    datasets: [
      {
        backgroundColor: 'rgba(0,0,0,0.2)',
        borderColor: 'rgba(0,0,0,1)',
        pointRadius: 0,
        data: [
          <?php
          $gh_contributors = (array) $stats_json->gh_contributors;
          sort($gh_contributors);
          $cumulative_count = 0;
          foreach($gh_contributors as $username => $timestamp){
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            $cumulative_count += 1;
            echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$cumulative_count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['gh_contribs'].options.title.text = 'nf-core GitHub code contributors over time';
  var ctx = document.getElementById('gh_contribs_plot').getContext('2d');
  charts['gh_contribs'] = new Chart(ctx, chartData['gh_contribs']);


  // Twitter followers chart
  chartData['twitter'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['twitter'].data = {
    datasets: [
      {
        backgroundColor: 'rgba(74, 161, 235, 0.2)',
        borderColor: 'rgba(74, 161, 235, 1)',
        pointRadius: 0,
        data: [
          <?php
          foreach($stats_json->twitter->followers_count as $timestamp => $count){
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['twitter'].options.title.text = '@nf_core twitter followers users over time';
  var ctx = document.getElementById('twitter_followers_plot').getContext('2d');
  charts['twitter'] = new Chart(ctx, chartData['twitter']);

  // Repo clones plot
  chartData['repo_clones'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['repo_clones'].data = {
    datasets: [
      {
        label: 'Clones per week',
        backgroundColor: 'rgba(83, 164, 81, 1.0)',
        borderColor: 'rgba(83, 164, 81, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-count',
        data: [
          <?php
          ksort($stats_total_allrepos['clones_count']);
          foreach($stats_total_allrepos['clones_count'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      },
      {
        label: 'Unique cloners per week',
        backgroundColor: 'rgba(33, 94, 190, 1.0)',
        borderColor: 'rgba(33, 94, 190, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-uniques',
        data: [
          <?php
          ksort($stats_total_allrepos['clones_uniques']);
          foreach($stats_total_allrepos['clones_uniques'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_clones'].options.title.text = 'nf-core git clones per week';
  chartData['repo_clones'].options.elements.line.borderWidth = 2;
  chartData['repo_clones'].options.scales.yAxes = [
    {
      id: 'y-axis-count',
      display: true,
      scaleLabel: {
        display: true,
        labelString: 'Clones per week',
        fontColor: 'rgba(83, 164, 81, 1.0)',
      },
      position: 'left',
      gridLines: { drawOnChartArea: false }
    },
    {
      id: 'y-axis-uniques',
      display: true,
      scaleLabel: {
        display: true,
        labelString: 'Unique cloners per week',
        fontColor: 'rgba(33, 94, 190, 1.0)',
      },
      position: 'right',
      gridLines: { drawOnChartArea: false }
    }
  ];
  var ctx = document.getElementById('repo_clones_plot').getContext('2d');
  charts['repo_clones'] = new Chart(ctx, chartData['repo_clones']);



  // Repo views plot
  chartData['repo_views'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['repo_views'].data = {
    datasets: [
      {
        label: 'Views per week',
        backgroundColor: 'rgba(83, 164, 81, 1.0)',
        borderColor: 'rgba(83, 164, 81, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-count',
        data: [
          <?php
          ksort($stats_total_allrepos['views_count']);
          foreach($stats_total_allrepos['views_count'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      },
      {
        label: 'Unique visitors per week',
        backgroundColor: 'rgba(33, 94, 190, 1.0)',
        borderColor: 'rgba(33, 94, 190, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-uniques',
        data: [
          <?php
          ksort($stats_total_allrepos['views_uniques']);
          foreach($stats_total_allrepos['views_uniques'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            // Skip zeros (anything before 2010)
            if($timestamp < 1262304000){
              continue;
            }
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_views'].options.title.text = 'nf-core repository web views per week';
  chartData['repo_views'].options.elements.line.borderWidth = 2;
  chartData['repo_views'].options.scales.yAxes = [
    {
      id: 'y-axis-count',
      display: true,
      scaleLabel: {
        display: true,
        labelString: 'Views per week',
        fontColor: 'rgba(83, 164, 81, 1.0)',
      },
      position: 'left',
      gridLines: { drawOnChartArea: false }
    },
    {
      id: 'y-axis-uniques',
      display: true,
      scaleLabel: {
        display: true,
        labelString: 'Unique visitors per week',
        fontColor: 'rgba(33, 94, 190, 1.0)',
      },
      position: 'right',
      gridLines: { drawOnChartArea: false }
    }
  ];
  var ctx = document.getElementById('repo_views_plot').getContext('2d');
  charts['repo_views'] = new Chart(ctx, chartData['repo_views']);

  // Make canvas2svg work with ChartJS
  // https://stackoverflow.com/a/52151467/713980
  function canvas2svgTweakLib(){
    C2S.prototype.getContext = function (contextId) {
      if (contextId=="2d" || contextId=="2D") { return this; }
      return null;
    }
    C2S.prototype.style = function () { return this.__canvas.style }
    C2S.prototype.getAttribute = function (name) { return this[name]; }
    C2S.prototype.addEventListener =  function(type, listener, eventListenerOptions) {
      console.log("canvas2svg.addEventListener() not implemented.")
    }
  }
  canvas2svgTweakLib();

  function exportChartJsSVG(target){
    // Turn off responiveness
    chartData[target].options.responsive = false;
    chartData[target].options.animation = false;
    // canvas2svg 'mock' context
    var svgContext = C2S(800,400);
    // new chart on 'mock' context fails:
    var mySvg = new Chart(svgContext, chartData[target]);
    // Failed to create chart: can't acquire context from the given item
    var svg = svgContext.getSerializedSvg(true);
    // Trigger browser download with SVG
    var blob = new Blob([svg], {
      type: "text/plain;charset=utf-8"
    });
    saveAs(blob, 'nf-core_'+target+'_plot.svg');
    // Turn responiveness back on again
    chartData[target].options.responsive = true;
    chartData[target].options.animation = true;
  }
  $('.dl_plot_svg').click(function(e){
    e.preventDefault();
    var target = $(this).data('target');
    exportChartJsSVG(target);
  });
  $('.reset_chart_zoom').click(function(e){
    e.preventDefault();
    var target = $(this).data('target');
    charts[target].resetZoom();
  });

});

</script>

<?php
$subfooter = '<p class="mb-0"><i class="far fa-clock"></i> Last updated: '.date('d-m-Y', $stats_json->updated).'</p>';

include('../includes/footer.php'); ?>
