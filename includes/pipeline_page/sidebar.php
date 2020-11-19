<?php
/////////////////////
// Sidebar for pipeline homepage with key stats
/////////////////////

// Get number of open issues and PRs
$issues_json_fn = dirname(dirname(dirname(__FILE__))).'/nfcore_issue_stats.json';
$issues_json = json_decode(file_get_contents($issues_json_fn), true);
$num_issues = 0;
foreach($issues_json['repos'][$pipeline->name]['issues'] as $issue){
  $num_issues += (int)($issue['state']=='open');
}

$num_prs = count($issues_json['repos'][$pipeline->name]['prs']);

// Get number of clones over time
$stats_json_fn = dirname(dirname(dirname(__FILE__))).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn), true);
$stats = $stats_json['pipelines'][$pipeline->name]['repo_metrics'][ $stats_json['updated'] ];
$clones_counts = $stats_json['pipelines'][$pipeline->name]['clones_count'];
$total_clones = 0;
$clones_since = false;
foreach($clones_counts as $datetime => $count){
  $total_clones += $count;
  if(!$clones_since) $clones_since = strtotime($datetime);
  $clones_since = min($clones_since, strtotime($datetime));
}

// Get contributor avatars
$contrib_avatars = [];
foreach($stats_json['pipelines'][$pipeline->name]['contributors'] as $contributor){
  $contributions = $contributor['total'];
  if ($contributions >1){
    $contributions .= " contributions";
  } else{
    $contributions .= " contribution";
  }
  $contrib_avatars[
    '<a href="'.$contributor['author']['html_url'].'" title="@'.$contributor['author']['login'].', '.$contributions.'" data-toggle="tooltip"><img src="'.$contributor['author']['avatar_url'].'"></a>'
  ] = $contributor['total'];
}
arsort($contrib_avatars);

// Last release and last commit
$last_release_time = 'N/A';
$release_cmd = ' -r '.$release;
if(count($pipeline->releases) > 0){
  $last_release_time = time_ago($pipeline->releases[0]->published_at);
}
$last_commit = time_ago($pipeline->updated_at);

?>

<div class="pipeline-sidebar">
  <div class="row border-bottom pb-2">
    <div class="col-12">
      <h6><i class="fas fa-terminal fa-xs"></i> command</h6>
      <div class="input-group input-group-sm pipeline-run-cmd">
        <input type="text" class="form-control input-sm code" id="pipeline-run-cmd-text" data-autoselect="" value="nextflow run <?php echo $pipeline->full_name; echo $release_cmd; ?> -profile test" aria-label="Copy run command" readonly="">
        <div class="input-group-append">
          <button class="btn btn-outline-secondary copy-txt" data-target="pipeline-run-cmd-text" data-toggle="tooltip" title="Copy to clipboard" type="button"><i class="fad fa-clipboard px-1"></i></button>
        </div>
      </div>
    </div>
  </div>

  <h6><i class="fas fa-arrow-down fa-xs"></i> <span id="clones_header">clones in last <?php echo time_ago($clones_since, false); ?></span></h6>
  <div class="row border-bottom">
    <div class="col-6">
      <p id="clones_count"><?php echo $total_clones; ?></p>
    </div>
    <div class="col-6" style="overflow: none;">
      <canvas id="clones_plot" height="70"></canvas>
    </div>
  </div>

  <div class="row border-bottom">
    <div class="col-6">
      <h6>stars</h6>
      <p><a href="<?php echo $pipeline->html_url;?>/stargazers"><?php echo $stats['stargazers_count']; ?></a></p>
    </div>
    <div class="col-6">
      <h6>watchers</h6>
      <p><a href="<?php echo $pipeline->html_url;?>/watchers"><?php echo $stats['subscribers_count']; ?></a></p>
    </div>
  </div>

  <div class="row border-bottom">
    <div class="col-6">
      <h6>last release</h6>
      <p><a href="/<?php echo $pipeline->name;?>/releases"><?php echo $last_release; ?></a></p>
    </div>
    <div class="col-6">
      <h6>last updated</h6>
      <p><?php echo $last_commit; ?></p>
    </div>
  </div>

  <div class="row border-bottom">
    <div class="col-6">
      <h6>open issues</h6>
      <p><a href="<?php echo $pipeline->html_url; ?>/issues"><?php echo $num_issues; ?></a></p>
    </div>
    <div class="col-6">
      <h6>pull requests</h6>
      <p><a href="<?php echo $pipeline->html_url; ?>/pulls"><?php echo $num_prs; ?></a></p>
    </div>
  </div>
  <div class="row border-bottom">
    <div class="col-12">
      <h6>collaborators</h6>
      <p class="contrib-avatars"><?php echo implode(array_keys($contrib_avatars)); ?></p>
    </div>
  </div>
  <div>
    <h6>get in touch</h6>
    <p><a class="btn btn-sm btn-outline-info" href="https://nfcore.slack.com/channels/<?php echo $pipeline->name; ?>"><i class="fab fa-slack mr-1"></i> Ask a question on Slack</a></p>
    <p><a class="btn btn-sm btn-outline-secondary" href="<?php echo $pipeline->html_url; ?>/issues"><i class="fab fa-github mr-1"></i> Open an issue on GitHub</a></p>
  </div>
</div>

<?php
// Collect this content into a variable to be inserted in to the very end of the HTML
ob_start();
?>

<div class="toast" id="pipeline_sidebar_cmd_copied" data-delay="5000" role="alert" aria-live="assertive" aria-atomic="true">
  <div class="toast-header">
    <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded mr-2" alt="">
    <strong class="mr-auto">Command copied!</strong>
    <button type="button" class="ml-2 mb-1 close" data-dismiss="toast" aria-label="Close">
      <span aria-hidden="true">&times;</span>
    </button>
  </div>
  <div class="toast-body">
    Paste this command into your terminal to run the pipeline with a small test dataset.
  </div>
</div>

<script type="text/javascript">
$(function(){
  // Plot hover vertical line
  // https://stackoverflow.com/a/45172506/713980
  Chart.defaults.LineWithLine = Chart.defaults.line;
  Chart.controllers.LineWithLine = Chart.controllers.line.extend({
     draw: function(ease) {
        Chart.controllers.line.prototype.draw.call(this, ease);

        if (this.chart.tooltip._active && this.chart.tooltip._active.length) {
           var activePoint = this.chart.tooltip._active[0],
               ctx = this.chart.ctx,
               x = activePoint.tooltipPosition().x,
               topY = this.chart.scales['y-axis-0'].top,
               bottomY = this.chart.scales['y-axis-0'].bottom;

           // draw line
           ctx.save();
           ctx.beginPath();
           ctx.moveTo(x, topY);
           ctx.lineTo(x, bottomY);
           ctx.lineWidth = 1;
           ctx.strokeStyle = '#999';
           ctx.stroke();
           ctx.restore();
        } else {
          $('#clones_header').text('clones in last <?php echo time_ago($clones_since, false); ?>');
          $('#clones_count').text('<?php echo $total_clones; ?>');
        }
     }
  });

  // Make the plot
  var ctx = document.getElementById('clones_plot').getContext('2d');
  new Chart(ctx, {
    data: {
      datasets: [
        {
          backgroundColor: 'rgba(84, 171, 106, 0.2)',
          borderColor: 'rgba(84, 171, 106, 1)',
          pointRadius: 0,
          pointHoverBorderColor: 'rgba(84, 171, 106, 0)', // transparent
          pointHoverBackgroundColor: 'rgba(84, 171, 106, 0)', // transparent
          data: [
            <?php
            $dates = [];
            foreach(array_keys($clones_counts) as $date){
              $dates[strtotime($date)] = $date;
            }
            ksort($dates);
            foreach($dates as $ts => $date){
              $count = $clones_counts[$date];
              echo '{ x: "'.date('Y-m-d', $ts).'", y: '.$count.' },'."\n\t\t\t";
            }
            ?>
            ]
          }
      ],
    },
    type: 'LineWithLine',
    options: {
      onClick: function(e){
        window.location.href = '/<?php echo $pipeline->name; ?>/stats';
      },
      elements: {
        point: {
          radius: 0,
          hitRadius: 3,
          hoverRadius: 3
        },
        line: {
          borderWidth: 2,
          tension: 0 // disables bezier curves
        }
      },
      scales: {
        xAxes: [{
          type: 'time',
          display: false
        }],
        yAxes: [{
          display: false
        }],
      },
      legend: {
        display: false
      },
      tooltips: {
        enabled: false,
        mode: 'x',
        intersect: false,
        custom: function(tooltipModel) {
          tooltipModel.opacity = 0
        },
        callbacks: {
          // Use the footer callback to display the sum of the items showing in the tooltip
          footer: function(tooltipItems, data) {
            $('#clones_header').text('clones - '+tooltipItems[0]['label']);
            $('#clones_count').text(tooltipItems[0]['value']);
          },
        }
      },
    }
  });
});
</script>

<?php
$end_of_html = ob_get_contents();
ob_end_clean();
