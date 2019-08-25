<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

$stats_json_fn = dirname(dirname(dirname(__FILE__))).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn), true);

ob_start();

?>

<h1><?php echo $pipeline->full_name; ?> statistics</h1>
<h2>Repository traffic</h2>
<div class="card mt-4">
  <div class="card-header">
    <span class="float-right small text-muted">
      <a href="#" data-target="repo_clones" class="dl_plot_svg text-muted"><i class="fas fa-download"></i> SVG</a>
      &nbsp;/&nbsp; <a href="#" data-target="repo_clones" class="reset_chart_zoom text-muted"><i class="fas fa-search-minus"></i> Reset zoom</a>
    </span>
    <a class="text-body" href="https://github.com/<?php echo $pipeline->full_name; ?>/graphs/traffic">Git clones</a>
  </div>
  <div class="card-body">
    <canvas id="repo_clones_plot" height="80"></canvas>
  </div>
  <div class="card-footer text-muted text-center small">
    <div class="row">
      <div class="col-6 border-right border-secondary">
        <span class="text-body lead"><?php echo round_nicely($stats_json['pipelines'][$pipeline->name]['clones_count_total']); ?></span>
        <br>Clones since <?php echo date('F Y', strtotime(array_keys($stats_json['pipelines'][$pipeline->name]['clones_count'])[0])); ?>
      </div>
      <div class="col-6">
        <span class="text-body lead"><?php echo round_nicely($stats_json['pipelines'][$pipeline->name]['clones_uniques_total']); ?></span>
        <br>Unique cloners since <?php echo date('F Y', strtotime(array_keys($stats_json['pipelines'][$pipeline->name]['clones_uniques'])[0])); ?>
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
    <a class="text-body" href="https://github.com/<?php echo $pipeline->full_name; ?>/graphs/traffic">Visitors</a>
  </div>
  <div class="card-body">
    <canvas id="repo_views_plot" height="80"></canvas>
  </div>
  <div class="card-footer text-muted text-center small">
    <div class="row align-items-center">
      <div class="col-6 border-right border-secondary">
        <span class="text-body lead"><?php echo round_nicely($stats_json['pipelines'][$pipeline->name]['views_count_total']); ?></span>
        <br>Views since <?php echo date('F Y', strtotime(array_keys($stats_json['pipelines'][$pipeline->name]['views_count'])[0])); ?>
      </div>
      <div class="col-6">
        <span class="text-body lead"><?php echo round_nicely($stats_json['pipelines'][$pipeline->name]['views_uniques_total']); ?></span>
        <br>Unique visitors since <?php echo date('F Y', strtotime(array_keys($stats_json['pipelines'][$pipeline->name]['views_uniques'])[0])); ?>
      </div>
    </div>
  </div>
</div>


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
          foreach($stats_json['pipelines'][$pipeline->name]['clones_count'] as $timestamp => $count){
            echo '{ x: "'.date('Y-m-d', strtotime($timestamp)).'", y: '.$count.' },'."\n\t\t\t";
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
          foreach($stats_json['pipelines'][$pipeline->name]['clones_uniques'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_clones'].options.title.text = 'nf-core git clones per day';
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
          foreach($stats_json['pipelines'][$pipeline->name]['views_count'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
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
          foreach($stats_json['pipelines'][$pipeline->name]['views_uniques'] as $timestamp => $count){
            $timestamp = strtotime($timestamp);
            echo '{ x: "'.date('Y-m-d', $timestamp).'", y: '.$count.' },'."\n\t\t\t";
          }
          ?>
        ]
      }
    ]
  };
  chartData['repo_views'].options.title.text = 'nf-core repository web views per day';
  chartData['repo_views'].options.scales.yAxes[0].scaleLabel.labelString = 'Views per day';
  chartData['repo_views'].options.scales.yAxes[1].scaleLabel.labelString = 'Unique visitors per day';
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

});

</script>

<?php

$content = ob_get_contents();
ob_end_clean();
