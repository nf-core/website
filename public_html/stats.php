<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$title = 'Community statistics';
$subtitle = 'Measuring activity across the nf-core community.';
include('../includes/header.php');

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

$pipeline_stats_json_fn = dirname(dirname(__FILE__)).'/nfcore_stats.json';
$pipeline_stats_json = json_decode(file_get_contents($pipeline_stats_json_fn));
$pipeline_stats = $pipeline_stats_json->pipelines;
# echo '<pre>'.print_r($pipeline_stats, true).'</pre>';

$stats_total = [
  'releases' => 0,
  'stargazers' => 0,
  'watchers' => 0,
  'forks' => 0,
  'clones_count_total' => 0,
  'clones_uniques_total' => 0,
  'views_count_total' => 0,
  'views_uniques_total' => 0,
  'unique_contributors' => []
];

$trows = [];
$missing_stats = [];
foreach($pipelines as $wf):
  $stats_total['releases'] += count($wf->releases);
  $stats_total['stargazers'] += $wf->stargazers_count;
  $stats_total['forks'] += $wf->forks_count;
  if(!isset($pipeline_stats->{$wf->name})){
    $missing_stats[] = $wf->full_name;
  } else {
    $stats_total['clones_count_total'] += $pipeline_stats->{$wf->name}->clones_count_total;
    $stats_total['clones_uniques_total'] += $pipeline_stats->{$wf->name}->clones_uniques_total;
    $stats_total['views_count_total'] += $pipeline_stats->{$wf->name}->views_count_total;
    $stats_total['views_uniques_total'] += $pipeline_stats->{$wf->name}->views_uniques_total;
    foreach($pipeline_stats->{$wf->name}->contributors as $contributor){
      $gh_username = $contributor->author->login;
      if(!isset($stats_total['unique_contributors'][$gh_username])){
        $stats_total['unique_contributors'][$gh_username] = 0;
      }
      $stats_total['unique_contributors'][$gh_username] += $contributor->total;
    }
  }
  $missing_stat = '<abbr title="New ppeline - stats not yet fetched" data-toggle="tooltip" class="bg-warning">&nbsp;?&nbsp;</abbr>';
  ob_start();
  ?>
  <tr>
    <td><?php echo '<a href="'.$wf->html_url.'" target="_blank">'.$wf->full_name.'</a>'; ?></td>
    <td><?php echo time_ago($wf->created_at, false); ?></td>
    <td class="text-right"><?php echo count($wf->releases); ?></td>
    <td class="text-right"><?php echo isset($pipeline_stats->{$wf->name}->num_contributors) ? $pipeline_stats->{$wf->name}->num_contributors : $missing_stat; ?></td>
    <td class="text-right"><?php echo $wf->stargazers_count; ?></td>
    <td class="text-right"><?php echo $wf->forks_count; ?></td>
    <td class="text-right"><?php echo isset($pipeline_stats->{$wf->name}->clones_count_total) ? $pipeline_stats->{$wf->name}->clones_count_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($pipeline_stats->{$wf->name}->clones_uniques_total) ? $pipeline_stats->{$wf->name}->clones_uniques_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($pipeline_stats->{$wf->name}->views_count_total) ? $pipeline_stats->{$wf->name}->views_count_total : $missing_stat; ?></td>
    <td class="text-right"><?php echo isset($pipeline_stats->{$wf->name}->views_uniques_total) ? $pipeline_stats->{$wf->name}->views_uniques_total : $missing_stat; ?></td>
  </tr>
<?php
$trows[] = ob_get_contents();
ob_end_clean();
endforeach;

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

<h1>Pipelines</h1>
<p>NB: Clones and repo view data collection started in July 2019.</p>


<table class="table table-sm small pipeline-stats-table">
  <thead class="thead-light">
    <tr>
      <th>Name</th>
      <th>Age</th>
      <th class="text-right">Releases</th>
      <th class="text-right"># Contributors</th>
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
      <th class="font-weight-light text-right"><?php echo $stats_total['releases']; ?></th>
      <th class="font-weight-light text-right"><?php echo count($stats_total['unique_contributors']); ?> unique</th>
      <th class="font-weight-light text-right"><?php echo $stats_total['stargazers']; ?></th>
      <th class="font-weight-light text-right"><?php echo $stats_total['forks']; ?></th>
      <th class="font-weight-light text-right"><?php echo $stats_total['clones_count_total']; ?></th>
      <th class="font-weight-light text-right"><?php echo $stats_total['clones_uniques_total']; ?></th>
      <th class="font-weight-light text-right"><?php echo $stats_total['views_count_total']; ?></th>
      <th class="font-weight-light text-right"><?php echo $stats_total['views_uniques_total']; ?></th>
    </tr>
  </thead>
  <tbody>
  <?php echo implode($trows); ?>
  </tbody>
  <tfoot class="thead-light">
    <tr>
      <th>Name</th>
      <th>Age</th>
      <th class="text-right">Releases</th>
      <th class="text-right"># Contributors</th>
      <th class="text-right">Stargazers</th>
      <th class="text-right">Forks</th>
      <th class="text-right">Clones</th>
      <th class="text-right">Unique cloners</th>
      <th class="text-right">Repo views</th>
      <th class="text-right">Unique repo visitors</th>
    </tr>
  </tfoot>
</table>

<script type="text/javascript">
$(function(){
  $('.toast').toast('show');
});
</script>

<?php
$subfooter = '<p class="mb-0"><i class="far fa-clock"></i> Last updated: '.date('h:i, d-m-Y', $pipeline_stats_json->updated).'</p>';

include('../includes/footer.php'); ?>
