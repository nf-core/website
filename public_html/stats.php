<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$title = 'Community statistics';
$subtitle = 'Measuring activity across the nf-core community.';
include('../includes/header.php');

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

$stats_json_fn = dirname(dirname(__FILE__)).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn));

// Convenience variables
$slack_users = $stats_json->slack->user_counts->{$stats_json->updated};
$twitter_users = $stats_json->twitter->followers_count->{$stats_json->updated};

# echo '<pre>'.print_r($stats, true).'</pre>';

// Run everything twice - keep pipelines and core repos seperate
foreach(['pipelines', 'core_repos'] as $repo_type):

$stats = $stats_json->{$repo_type};
$stats_total[$repo_type] = [
  'releases' => 0,
  'stargazers' => 0,
  'watchers' => 0,
  'forks' => 0,
  'clones_count_total' => 0,
  'clones_uniques_total' => 0,
  'views_count_total' => 0,
  'views_uniques_total' => 0,
  'unique_contributors' => [],
  'total_commits' => 0
];

$trows[$repo_type] = [];
$missing_stats = [];
foreach($stats as $repo_name => $repo):
  $metrics = $repo->repo_metrics->{$stats_json->updated};
  $releases = 0;
  foreach($pipelines as $wf){
    if($wf->name == $repo_name){
      $releases = count($wf->releases);
    }
  }
  $stats_total[$repo_type]['releases'] += $releases;
  $stats_total[$repo_type]['stargazers'] += $metrics->stargazers_count;
  $stats_total[$repo_type]['forks'] += $metrics->forks_count;
  $total_commits = 0;
  if(!isset($repo)){
    $missing_stats[] = $metrics->full_name;
  } else {
    $stats_total[$repo_type]['clones_count_total'] += $repo->clones_count_total;
    $stats_total[$repo_type]['clones_uniques_total'] += $repo->clones_uniques_total;
    $stats_total[$repo_type]['views_count_total'] += $repo->views_count_total;
    $stats_total[$repo_type]['views_uniques_total'] += $repo->views_uniques_total;
    foreach($repo->contributors as $contributor){
      $gh_username = $contributor->author->login;
      $stats_total[$repo_type]['unique_contributors'][$gh_username] = 0;
      $stats_total['total']['unique_contributors'][$gh_username] = 0;
      $stats_total[$repo_type]['total_commits'] += $contributor->total;
      $total_commits += $contributor->total;
    }
  }
  $missing_stat = '<abbr title="New ppeline - stats not yet fetched" data-toggle="tooltip" class="bg-warning">&nbsp;?&nbsp;</abbr>';
  ob_start();
  ?>
  <tr>
    <td><?php echo '<a href="'.$metrics->html_url.'" target="_blank">'.$metrics->full_name.'</a>'; ?></td>
    <td><?php echo time_ago($metrics->created_at, false); ?></td>
    <?php if($repo_type == 'pipelines'): ?><td class="text-right"><?php echo $releases; ?></td><?php endif; ?>
    <td class="text-right"><?php echo isset($repo->num_contributors) ? $repo->num_contributors : $missing_stat; ?></td>
    <td class="text-right"><?php echo $total_commits; ?></td>
    <td class="text-right"><?php echo $metrics->stargazers_count; ?></td>
    <td class="text-right"><?php echo $metrics->forks_count; ?></td>
    <td class="text-right"><?php echo isset($repo->clones_count_total) ? $repo->clones_count_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($repo->clones_uniques_total) ? $repo->clones_uniques_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($repo->views_count_total) ? $repo->views_count_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($repo->views_uniques_total) ? $repo->views_uniques_total : $missing_stat; ?></td>
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

<h1>nf-core in numbers</h1>
<p class="text-info small">
  <i class="far fa-hand-point-right"></i>
  Click a number to see how the community has grown over time
</p>

<div class="card-group text-center stats_keynumbers">
  <div class="card bg-light" data-toggle="collapse" data-target="#slack_chart">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $slack_users->total; ?></p>
      <p class="card-text text-muted">Slack users</p>
    </div>
    <div class="bg-icon" style="color: rgba(89, 37, 101, 0.1);"><i class="fab fa-slack"></i></div>
  </div>
  <div class="card bg-light" data-toggle="collapse" data-target="#gh_orgmembers_chart">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $stats_json->gh_org_members->{$stats_json->updated}; ?></p>
      <p class="card-text text-muted">GitHub organisation members</p>
    </div>
    <div class="bg-icon"><i class="fab fa-github"></i></div>
  </div>
  <div class="card bg-light" data-toggle="collapse" data-target="#gh_contribs_chart">
    <div class="card-body">
      <p class="card-text display-4"><?php echo count($stats_total['total']['unique_contributors']); ?></p>
      <p class="card-text text-muted">GitHub contributors</p>
    </div>
    <div class="bg-icon"><i class="fas fa-code-branch"></i></div>
  </div>
  <div class="card bg-light" data-toggle="collapse" data-target="#twitter_chart">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $twitter_users; ?></p>
      <p class="card-text text-muted">Twitter followers</p>
    </div>
    <div class="bg-icon" style="color: rgba(74, 161, 235, 0.2);"><i class="fab fa-twitter"></i></div>
  </div>
</div>
<div id="stats_keynumbers_chart_wrapper">
  <div class="card bg-light collapse stats_keynumbers_chart" id="slack_chart" data-parent="#stats_keynumbers_chart_wrapper">
    <div class="card-body">
      <canvas id="slack_users_plot" width="400" height="130"></canvas>
      <p class="card-text small text-muted text-center"><i class="fas fa-info-circle"></i> Slack considers users to be active when they haven't used slack for the previous 14 days.</p>
    </div>
  </div>
  <div class="card bg-light collapse stats_keynumbers_chart" id="gh_orgmembers_chart" data-parent="#stats_keynumbers_chart_wrapper">
    <div class="card-body">
      GitHub org members chart
    </div>
  </div>
  <div class="card bg-light collapse stats_keynumbers_chart" id="gh_contribs_chart" data-parent="#stats_keynumbers_chart_wrapper">
    <div class="card-body">
      GitHub contribs chart
    </div>
  </div>
  <div class="card bg-light collapse stats_keynumbers_chart" id="twitter_chart" data-parent="#stats_keynumbers_chart_wrapper">
    <div class="card-body">
      <canvas id="twitter_followers_plot" width="400" height="130"></canvas>
    </div>
  </div>
</div>




<?php
// The pipeline and core repo tables are the same
foreach(['pipelines', 'core_repos'] as $repo_type):

if(count($missing_stats)){ echo '<div class="toasts-container">'; }
foreach($missing_stats as $missing_stat): ?>
<div class="toast" role="alert" aria-live="assertive" aria-atomic="true" data-delay="5000">
  <div class="toast-header">
    <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded mr-2" alt="nf-core logo">
    <strong class="mr-auto">nf-core</strong>
    <button type="button" class="ml-2 mb-1 close" data-dismiss="toast" aria-label="Close">
      <span aria-hidden="true">&times;</span>
    </button>
  </div>
  <div class="toast-body">
    Could not fetch statistics for <code><?php echo $missing_stat; ?></code>
  </div>
</div>
<?php
endforeach;
if(count($missing_stats)){ echo '</div>'; }
?>

<h1><?php echo ucfirst(str_replace('_', ' ', $repo_type)); ?></h1>

<div class="card mb-3">
  <div class="card-header">
      <button class="btn btn-link p-0 pb-1 collapsed text-muted" type="button" data-toggle="collapse" data-target="#caveats_<?php echo $repo_type; ?>">
        Read about how these numbers are collected and what caveats should be considered
      </button>
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
            <li>Traditional HPC centres may share workflow installations, so only have one clone for many users / pipeline runs</li>
            <li>Cloud users will typically spin up a new instance and clone the workflow every time that they run a pipeline.</li>
          </ul>
        </li>
        <li>Clone counts and repositoriy views are only available for two weeks - longer term data collection for nf-core repos started in July 2019. This is when we started counting the totals.</li>
        <li>Metrics are fetched using the GitHub API only once per week (last checked <?php echo date('d-m-Y', $stats_json->updated); ?>).</li>
      </ul>
    </div>
  </div>
</div>


<table class="table table-sm small pipeline-stats-table">
  <thead class="thead-light">
    <tr>
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

<?php endforeach;

?>

<script type="text/javascript">
$(function(){
  $('.toast').toast('show');

  // Slack users chart
  var ctx = document.getElementById('slack_users_plot').getContext('2d');
  var slack_users_plot = new Chart(ctx, {
    type: 'line',
    data: {
      datasets: [
        {
          label: 'Inactive',
          backgroundColor: 'rgba(150,150,150, 0.2)',
          borderColor: 'rgba(150,150,150, 1)',
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
          data: [
            <?php
            foreach($stats_json->slack->user_counts as $timestamp => $users){
              echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$users->active.' },'."\n\t\t\t";
            }
            ?>
          ]
        }
      ]
    },
    options: {
      title: {
        display: true,
        text: 'nf-core Slack users over time'
      },
      elements: {
        line: {
          borderWidth: 1,
          tension: 0 // disables bezier curves
        }
      },
      scales: {
        xAxes: [{
          type: 'time'
        }],
        yAxes: [{
          stacked: true
        }]
      },
      legend: {
        position: 'bottom',
        labels: {
          lineWidth: 1
        }
      },
      tooltips: {
        mode: 'x'
      },
    }
  });


  // Twitter followers chart
  var ctx = document.getElementById('twitter_followers_plot').getContext('2d');
  var twitter_followers_plot = new Chart(ctx, {
    type: 'line',
    data: {
      datasets: [
        {
          backgroundColor: 'rgba(74, 161, 235, 0.2)',
          borderColor: 'rgba(74, 161, 235, 1)',
          data: [
            <?php
            foreach($stats_json->twitter->followers_count as $timestamp => $count){
              echo '{ x: "'.date('Y-m-d H:i:s', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
            }
            ?>
          ]
        }
      ]
    },
    options: {
      title: {
        display: true,
        text: '@nf_core twitter followers users over time'
      },
      elements: {
        line: {
          borderWidth: 1,
          tension: 0 // disables bezier curves
        }
      },
      scales: {
        xAxes: [{
          type: 'time'
        }],
        yAxes: [{
          stacked: true
        }]
      },
      legend: {
        display: false
      },
      tooltips: {
        mode: 'x'
      },
    }
  });
});

</script>

<?php
$subfooter = '<p class="mb-0"><i class="far fa-clock"></i> Last updated: '.date('d-m-Y', $stats_json->updated).'</p>';

include('../includes/footer.php'); ?>
