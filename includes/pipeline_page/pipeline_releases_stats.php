<?php
// Build the HTML for a pipeline documentation page.
// Imported by pipeline_page/_index.php

// Load pipeline stats
$stats_json_fn = dirname(dirname(dirname(__FILE__))).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn), true);
$stats = $stats_json['pipelines'][$pipeline->name];
$metrics = $stats['repo_metrics'][ $stats_json['updated'] ];

// Load full contributor stats for this pipeline
$contrib_json_fn = dirname(dirname(dirname(__FILE__))).'/contributor_stats/'.$pipeline->name.'.json';
$contrib_json = json_decode(file_get_contents($contrib_json_fn), true);

ob_start();
echo '<h1 class="mt-0">Version history</h1>';

$first = true;
foreach($pipeline->releases as $releases){ ?>

<div class="row">
  <div class="col-auto">
    <a href="#download-<?php echo $releases->tag_sha; ?>" class="text-body" data-toggle="collapse">
      <samp><?php echo $releases->tag_name; ?></samp>
    </a>
  </div>
  <div class="col">
  </div>
  <div class="col-auto">
    <a href="#download-<?php echo $releases->tag_sha; ?>" class="text-body" data-toggle="collapse"><small class="text-muted"><?php echo time_ago($releases->published_at); ?></small></a>
    <button class="btn btn-sm btn-link text-body" type="button" data-toggle="collapse" data-target="#download-<?php echo $releases->tag_sha; ?>">
      <i class="fas fa-caret-<?php echo $first ? 'down' : 'left'; ?>"></i>
  </button>
  </div>
</div>
<div class="collapse <?php if($first){ echo 'show'; } ?>" id="download-<?php echo $releases->tag_sha; ?>">
  <div class="row pb-2">
    <div class="col-sm-6 small">
      Released <?php echo date('j M Y', strtotime($releases->published_at)); ?> &mdash;
      <code><?php echo substr($releases->tag_sha, 0, 7); ?></code>
    </div>
    <div class="col-sm-6 text-right">
      <a href="<?php echo $releases->zipball_url; ?>" class="btn btn-sm btn-outline-success">Download .zip</a>
      <a href="<?php echo $releases->tarball_url; ?>" class="btn btn-sm btn-outline-success">Download .tar.gz</a>
      <a href="<?php echo $releases->html_url; ?>" class="btn btn-sm btn-success"><i class="fab fa-github"></i> View release</a>
    </div>
  </div>
</div>
<hr class="m-0">



<?php $first = false;
}
?>

<script type="text/javascript">
$(function(){
  $('.collapse').on('show.bs.collapse', function () {
    var target = $(this).attr('id');
    $("button[data-target='#" + target + "']").html('<i class="fas fa-caret-down"></i>');
  });
  $('.collapse').on('hide.bs.collapse', function () {
    var target = $(this).attr('id');
    $("button[data-target='#" + target + "']").html('<i class="fas fa-caret-left"></i>');
  });
});
</script>

<h1 id="stats" class="mt-5"><a href="#stats" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a> Pipeline Statistics</h1>
<div class="card-group text-center stats_keynumbers mt-5 pb-5">
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $metrics['stargazers_count']; ?></p>
      <p class="card-text text-muted">Stars</p>
    </div>
    <div class="bg-icon"><i class="far fa-star"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#github_prs" class="text-body text-decoration-none stretched-link"><?php echo $metrics['network_forks_count']; ?></a></p>
      <p class="card-text text-muted">Forks</p>
    </div>
    <div class="bg-icon"><i class="fas fa-code-branch"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><?php echo $stats['commits']; ?></p>
      <p class="card-text text-muted">Commits</p>
    </div>
    <div class="bg-icon"><i class="far fa-file-code"></i></div>
  </div>
  <div class="card bg-light">
    <div class="card-body">
      <p class="card-text display-4"><a href="#github_issues" class="text-body text-decoration-none stretched-link"><?php echo $stats['num_contributors']; ?></a></p>
      <p class="card-text text-muted">Code contributors</p>
    </div>
    <div class="bg-icon"><i class="fas fa-user"></i></div>
  </div>
</div>

<h2 id="traffic"><a href="#traffic" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Repository traffic</h2>
<div class="card mt-4">
  <div class="card-header">
    <span class="float-right small text-muted">
      <a href="#" data-target="repo_clones" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
      &nbsp;/&nbsp; <a href="#" data-target="repo_clones" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
    </span>
    <a class="text-body" href="https://github.com/<?php echo $pipeline->full_name; ?>/graphs/traffic">Git clones</a>
  </div>
  <div class="card-body">
    <div style="height: 250px;"><canvas id="repo_clones_plot" height="80"></canvas></div>
  </div>
  <div class="card-footer text-muted text-center small">
    <div class="row">
      <div class="col-6 border-right border-secondary">
        <span class="text-body lead"><?php echo round_nicely($stats['clones_count_total']); ?></span>
        <br>Clones since <?php echo date('F Y', strtotime(array_keys($stats['clones_count'])[0])); ?>
      </div>
      <div class="col-6">
        <span class="text-body lead"><?php echo round_nicely($stats['clones_uniques_total']); ?></span>
        <br>Unique cloners since <?php echo date('F Y', strtotime(array_keys($stats['clones_uniques'])[0])); ?>
      </div>
    </div>
  </div>
</div>

<div class="card mt-4 mb-5">
  <div class="card-header">
    <span class="float-right small text-muted">
      <a href="#" data-target="repo_views" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
      &nbsp;/&nbsp; <a href="#" data-target="repo_views" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
    </span>
    <a class="text-body" href="https://github.com/<?php echo $pipeline->full_name; ?>/graphs/traffic">Visitors</a>
  </div>
  <div class="card-body">
    <div style="height: 250px;"><canvas id="repo_views_plot" height="80"></canvas></div>
  </div>
  <div class="card-footer text-muted text-center small">
    <div class="row align-items-center">
      <div class="col-6 border-right border-secondary">
        <span class="text-body lead"><?php echo round_nicely($stats['views_count_total']); ?></span>
        <br>Views since <?php echo date('F Y', strtotime(array_keys($stats['views_count'])[0])); ?>
      </div>
      <div class="col-6">
        <span class="text-body lead"><?php echo round_nicely($stats['views_uniques_total']); ?></span>
        <br>Unique visitors since <?php echo date('F Y', strtotime(array_keys($stats['views_uniques'])[0])); ?>
      </div>
    </div>
  </div>
</div>


<h2 id="contributors"><a href="#contributors" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Contributors</h2>

<div style="height: 250px;"><canvas id="contributors_plot" height="80"></canvas></div>

<div class="alert alert-info small p-2 mt-3 mb-3" role="alert">
  <i class="far fa-hand-point-right pl-2 pr-2"></i>
  Hover over the plot or an author's avatar to highlight commits
</div>
<p class="contrib-avatars">
<?php
$contrib_avatars = [];
foreach($contrib_json as $contrib){
  $contrib_avatars[
    '<a class="d-inline-block" href="https://github.com/'.$pipeline->full_name.'/graphs/contributors" data-author="'.$contrib['author']['login'].'" data-toggle="tooltip" title="@'.$contrib['author']['login'].'"><img src="'.$contrib['author']['avatar_url'].'"></a>'
  ] = $contrib['total'];
}
arsort($contrib_avatars);
echo implode(array_keys($contrib_avatars));
?>
</p>

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
          borderWidth: 2,
          tension: 0 // disables bezier curves
        }
      },
      maintainAspectRatio: false,
      scales: {
        xAxes: [{
          type: 'time',
          time: { minUnit: 'day' }
        }],
        yAxes: [
          {
            id: 'y-axis-count',
            display: true,
            scaleLabel: {
              display: true,
              labelString: 'Clones per day',
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
              labelString: 'Unique cloners per day',
              fontColor: 'rgba(33, 94, 190, 1.0)',
            },
            position: 'right',
            gridLines: { drawOnChartArea: false }
          }
        ]
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

  // Repo clones plot
  chartData['repo_clones'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['repo_clones'].data = {
    datasets: [
      {
        label: 'Clones per day',
        backgroundColor: 'rgba(83, 164, 81, 1.0)',
        borderColor: 'rgba(83, 164, 81, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-count',
        data: [
          <?php
          $dates = [];
          foreach(array_keys($stats['clones_count']) as $date){
            $dates[strtotime($date)] = $date;
          }
          ksort($dates);
          foreach($dates as $ts => $date){
            $count = $stats['clones_count'][$date];
            echo '{ x: "'.date('Y-m-d', $ts).'", y: '.$count.' },'."\n\t\t\t";
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
          $dates = [];
          foreach(array_keys($stats['clones_uniques']) as $date){
            $dates[strtotime($date)] = $date;
          }
          ksort($dates);
          foreach($dates as $ts => $date){
            $count = $stats['clones_uniques'][$date];
            echo '{ x: "'.date('Y-m-d', $ts).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_clones'].options.title.text = '<?php echo $pipeline->full_name; ?> git clones per day';
  chartData['repo_clones'].options.scales.yAxes[0].scaleLabel.labelString = 'Clones per day';
  chartData['repo_clones'].options.scales.yAxes[1].scaleLabel.labelString = 'Unique cloners per day';
  var ctx = document.getElementById('repo_clones_plot').getContext('2d');
  charts['repo_clones'] = new Chart(ctx, chartData['repo_clones']);



  // Repo views plot
  chartData['repo_views'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['repo_views'].data = {
    datasets: [
      {
        label: 'Views per day',
        backgroundColor: 'rgba(83, 164, 81, 1.0)',
        borderColor: 'rgba(83, 164, 81, 1.0)',
        fill: false,
        pointRadius: 0,
        yAxisID: 'y-axis-count',
        data: [
          <?php
          $dates = [];
          foreach(array_keys($stats['views_count']) as $date){
            $dates[strtotime($date)] = $date;
          }
          ksort($dates);
          foreach($dates as $ts => $date){
            $count = $stats['views_count'][$date];
            echo '{ x: "'.date('Y-m-d', $ts).'", y: '.$count.' },'."\n\t\t\t";
          }
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
          $dates = [];
          foreach(array_keys($stats['views_uniques']) as $date){
            $dates[strtotime($date)] = $date;
          }
          ksort($dates);
          foreach($dates as $ts => $date){
            $count = $stats['views_uniques'][$date];
            echo '{ x: "'.date('Y-m-d', $ts).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_views'].options.title.text = '<?php echo $pipeline->full_name; ?> repository web views per day';
  chartData['repo_views'].options.scales.yAxes[0].scaleLabel.labelString = 'Views per day';
  chartData['repo_views'].options.scales.yAxes[1].scaleLabel.labelString = 'Unique visitors per day';
  var ctx = document.getElementById('repo_views_plot').getContext('2d');
  charts['repo_views'] = new Chart(ctx, chartData['repo_views']);




  // Contributors plot
  chartData['contributors'] = JSON.parse(JSON.stringify(chartjs_base));
  chartData['contributors'].data = {
    datasets: [
      <?php $first = true;
      foreach($contrib_json as $contrib){ ?>
        {
          label: '<?php echo $contrib['author']['login']; ?>',
          backgroundColor: 'rgba(200, 200, 200, 0.5)',
          pointRadius: 0,
          <?php if($first) { echo "fill: 'origin',"; $first = false; } ?>
          data: [
            <?php
            foreach($contrib['weeks'] as $week){
              echo '{ x: "'.date('Y-m-d', $week['w']).'", y: '.$week['c'].' },'."\n\t\t\t";
            }
            ?>
          ]
        },
      <?php } ?>
    ]
  };
  chartData['contributors'].options.title.text = '<?php echo $pipeline->full_name; ?> commits per week';
  chartData['contributors'].options.scales.xAxes[0].gridLines = { display: false };
  chartData['contributors'].options.scales.yAxes = [{
    stacked: true,
    scaleLabel: {
      display: true,
      labelString: '# Commits per week'
    }
  }];
  chartData['contributors'].options.elements.line.fill = '-1'; // by default, fill lines to the previous dataset
  chartData['contributors'].options.elements.line.borderWidth = 0;
  chartData['contributors'].options.tooltips = {
    enabled: false,
    mode: 'nearest',
    intersect: false,
    custom: function(tooltipModel) {
      tooltipModel.opacity = 0
    },
    callbacks: {
      // Use the footer callback to display the sum of the items showing in the tooltip
      footer: function(tooltipItems, data) {
        var dsidx = tooltipItems[0].datasetIndex;
        var author = charts['contributors'].data.datasets[dsidx].label;
        // Highlight avatar
        $('.contrib-avatars a').tooltip('hide');
        $('.contrib-avatars a img').css({
          'filter': 'grayscale(100%)',
          'opacity': '0.2',
        });
        $(".contrib-avatars a[data-author='"+author+"']").tooltip('show');
        $(".contrib-avatars a[data-author='"+author+"'] img").css({
          'filter': 'grayscale(0%)',
          'opacity': '1',
        });
        // Higlight plot series
        $.each(charts['contributors'].data.datasets, function( idx, dataset ) {
          if(idx == dsidx){
            dataset.backgroundColor = '#22ae63';
          } else {
            dataset.backgroundColor = 'rgba(200, 200, 200, 0.5)';
          }
        });
        charts['contributors'].update(0);
      },
    }
  };
  var ctx = document.getElementById('contributors_plot').getContext('2d');
  charts['contributors'] = new Chart(ctx, chartData['contributors']);


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
    var svgContext = C2S(canvas_width,canvas_height);
    // new chart on 'mock' context fails:
    var mySvg = new Chart(svgContext, chartData[target]);
    // Failed to create chart: can't acquire context from the given item
    var svg = svgContext.getSerializedSvg(true);
    // Trigger browser download with SVG
    var blob = new Blob([svg], {
      type: "text/plain;charset=utf-8"
    });
    saveAs(blob, 'nf-core_<?php echo $pipeline->name; ?>_'+target+'_plot.svg');

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


  $('.contrib-avatars a').hover(function(){
    // Highlight avatar
    $('.contrib-avatars a img').css({
      'filter': 'grayscale(100%)',
      'opacity': '0.2',
    });
    $(this).find('img').css({
      'filter': 'grayscale(0%)',
      'opacity': '1',
    });
    // Higlight plot series
    var author = $(this).data('author');
    $.each(charts['contributors'].data.datasets, function( idx, dataset ) {
      if(dataset.label == author){
        dataset.backgroundColor = '#22ae63';
      } else {
        dataset.backgroundColor = 'rgba(200, 200, 200, 0.5)';
      }
    });
    charts['contributors'].update(0);
  });
  $('.contrib-avatars').mouseout(function(){
    $('.contrib-avatars a img').css({
      'filter': 'grayscale(0%)',
      'opacity': '1',
    });
    $.each(charts['contributors'].data.datasets, function( idx, dataset ) {
      dataset.backgroundColor = 'rgba(200, 200, 200, 0.5)';
    });
    charts['contributors'].update(0);
  });

});

</script>

<?php

$content = '<div class="pipeline-stats">'.ob_get_contents().'</div>';
ob_end_clean();
