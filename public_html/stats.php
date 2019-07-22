<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$title = 'nf-core statistics';
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
foreach($pipelines as $wf):
  $stats_total['releases'] += count($wf->releases);
  $stats_total['stargazers'] += $wf->stargazers_count;
  $stats_total['forks'] += $wf->forks_count;
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
  ob_start();
  ?>
  <tr>
    <td><?php echo '<a href="'.$wf->html_url.'" target="_blank">'.$wf->full_name.'</a>'; ?></td>
    <td><?php echo time_ago($wf->created_at, false); ?></td>
    <td class="text-right"><?php echo count($wf->releases); ?></td>
    <td class="text-right"><?php echo $pipeline_stats->{$wf->name}->num_contributors; ?></td>
    <td class="text-right"><?php echo $wf->stargazers_count; ?></td>
    <td class="text-right"><?php echo $wf->forks_count; ?></td>
    <td class="text-right"><?php echo $pipeline_stats->{$wf->name}->clones_count_total; ?></td>
    <td class="text-right"><?php echo $pipeline_stats->{$wf->name}->clones_uniques_total; ?></td>
    <td class="text-right"><?php echo $pipeline_stats->{$wf->name}->views_count_total; ?></td>
    <td class="text-right"><?php echo $pipeline_stats->{$wf->name}->views_uniques_total; ?></td>
  </tr>
<?php
$trows[] = ob_get_contents();
ob_end_clean();
endforeach;
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

<?php include('../includes/footer.php'); ?>
