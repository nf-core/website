<?php
/////////////////////
// Sidebar for pipeline homepage with key stats
/////////////////////

// // Get number of clones over time
// $stats_json_fn = dirname(dirname(dirname(__FILE__))) . '/nfcore_stats.json';
// $stats_json = json_decode(file_get_contents($stats_json_fn), true);
// $stats = $stats_json['pipelines'][$pipeline->name]['repo_metrics'][$stats_json['updated']];

$total_clones = array_sum(array_column($traffic_stats, 'clones'));
$clones_since = end($traffic_stats)['timestamp'];

// Get contributor avatars
$contrib_avatars = [];
foreach (array_unique(array_column($contributor_stats, 'author', 'avatar_url')) as $contributor_url => $contributor) {
    // print_r($contributor."\n");
    // get all contributions for this contributor
    $contributions = array_filter($contributor_stats, function ($c) use ($contributor) {
        return $c['author'] === $contributor;
    });
    // count the contributions
    $total_commits = array_sum(array_column($contributions, 'week_commits'));
    if ($total_commits > 1) {
        $total_commits .= ' contributions';
    } else {
        $total_commits .= ' contribution';
    }
    $contrib_avatars[
        '<a href="https://github.com/' .
            $contributor .
            '" title="@' .
            $contributor .
            ', ' .
            $total_commits .
            '" data-bs-toggle="tooltip"><img src="' .
            $contributor_url .
            '"></a>'
    ] = $total_commits;
}
arsort($contrib_avatars);

// Last release and last commit
$last_release_time = 'N/A';
$release_cmd = ' -r ' . $release;
if (!is_null($pipeline_metrics['last_release_date'])) {
    $last_release_time = time_ago($pipeline_metrics['last_release_date']);
}
$last_commit = time_ago($pipeline_metrics['gh_updated_at']);

// filter events for embed_at key and look if it is set to the this pipeline
$pipeline_name = explode('/', $_GET['path'])[0];
$embed_video = array_filter($events, function ($event) {
    global $pipeline_name;
    if (array_key_exists('embed_at', $event) && $pipeline_name === $event['embed_at']) {
        return $event;
    }
});
$embed_video = array_values($embed_video)[0];
?>

<div class="pipeline-sidebar">
    <div class="row pb-2 mb-md-4">
        <div class="col-12">
            <h6><i class="fas fa-terminal fa-xs"></i> Run with</h6>
            <ul class="nav nav-tabs border-bottom-0" id="myTab" role="tablist">
                <li class="nav-item" role="presentation">
                    <button class="nav-link text-muted active" id="nfcore-tab" data-bs-toggle="tab" data-bs-target="#pipeline-nfcore-run-cmd" type="button" role="tab" aria-controls="nfcore" aria-selected="true">nf-core</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link text-muted" id="nf-tab" data-bs-toggle="tab" data-bs-target="#nf" type="button" role="tab" aria-controls="nf" aria-selected="false">Nextflow</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link text-muted" id="tw-tab" data-bs-toggle="tab" data-bs-target="#tw" type="button" role="tab" aria-controls="tw" aria-selected="false">Tower</button>
                </li>
            </ul>
            <div class="tab-content mt-2">
                <div class="tab-pane show active" id="pipeline-nfcore-run-cmd" role="tabpanel" aria-labelledby="nfcore-tab">
                    <div class=" input-group input-group-sm pipeline-run-cmd">
                        <input type="text" class="form-control input-sm code rounded-0" id="pipeline-nfcore-run-cmd-text"  value="nf-core launch <?php echo $pipeline->full_name .
                            $release_cmd; ?>" aria-label="Copy run command" readonly="">
                        <button class="btn btn-outline-secondary copy-txt rounded-0" data-bs-target="pipeline-nfcore-run-cmd-text" data-bs-toggle="tooltip" data-bs-placement="left" title="Copy to clipboard" type="button"><i class="fas fa-clipboard px-1"></i></button>
                    </div>
                </div>
                <div class="tab-pane" id="nf" role="tabpanel" aria-labelledby="nf-tab">
                    <div class=" input-group input-group-sm pipeline-run-cmd">
                        <input type="text" class="form-control input-sm code  rounded-0" id="pipeline-nf-run-cmd-text"  value="nextflow run <?php echo $pipeline->full_name .
                            $release_cmd; ?> -profile test --outdir <OUTDIR>" aria-label="Copy run command" readonly="">
                        <button class="btn btn-outline-secondary copy-txt rounded-0" data-bs-target="pipeline-nf-run-cmd-text" data-bs-toggle="tooltip" data-bs-placement="left" title="Copy to clipboard" type="button"><i class="fas fa-clipboard px-1"></i></button>
                    </div>
                </div>
                <div class="tab-pane" id="tw" role="tabpanel" aria-labelledby="tw-tab">
                    <div class=" input-group input-group-sm pipeline-run-cmd">
                        <input type="text" class="form-control input-sm code  rounded-0" id="pipeline-tw-run-cmd-text" data-autoselect="" value="tw launch https://nf-co.re/<?php echo $pipeline->name .
                            $release_cmd; ?>" aria-label="Copy run command" readonly="">
                            <button class="btn btn-outline-secondary copy-txt rounded-0" data-bs-target="pipeline-tw-run-cmd-text" data-bs-toggle="tooltip" data-bs-placement="left" title="Copy to clipboard" type="button"><i class="fas fa-clipboard px-1"></i></button>
                </div><p class="text-muted">Read how to configure the Tower CLI <u><a href='https://github.com/seqeralabs/tower-cli/#2-configuration' target="_blank">here</a></u>.</p></div>
            </div>
        </div>
    </div>
    <?php if (isset($embed_video)): ?>
        <div class="row border-bottom">
            <div class="col-12">
                <h6><i class="fab fa-youtube fa-xs"></i> video introduction</h6>
                <div class="ratio ratio-16x9 mb-2">
                    <?php
                    //taken from  https://stackoverflow.com/a/8260383
                    preg_match(
                        '/^.*((youtu.be\/)|(v\/)|(\/u\/\w\/)|(embed\/)|(watch\?))\??v?=?([^#&?]*).*/',
                        $embed_video['youtube_embed'],
                        $matches,
                    );
                    $youtube_id = $matches[7];
                    ?>
                    <!-- Using the lazy loading trick from https://css-tricks.com/lazy-load-embedded-youtube-videos/-->
                    <iframe src="https://www.youtube.com/embed/<?php echo $youtube_id; ?>?autoplay=1" srcdoc="<style>*{padding:0;margin:0;overflow:hidden}html,body{height:100%}img,span{position:absolute;width:100%;top:0;bottom:0;margin:auto}span{height:1.5em;text-align:center;font:48px/1.5 sans-serif;color:white;text-shadow:0 0 0.5em black}</style>
                        <a href=https://www.youtube.com/embed/<?php echo $youtube_id; ?>?autoplay=1><img src=https://img.youtube.com/vi/<?php echo $youtube_id; ?>/hqdefault.jpg alt=''>
                            <span>▶</span>
                        </a>" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen title="">
                    </iframe>
                </div>
            </div>
        </div>
    <?php endif; ?>
    <h6><i class="fas fa-arrow-down fa-xs"></i> <span id="clones_header">clones in last
            <?php echo time_ago($clones_since, false); ?></span></h6>
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
            <p><a href="<?php echo $pipeline->html_url; ?>/stargazers"><?php echo $pipeline_metrics[
    'stargazers_count'
]; ?></a></p>
        </div>
        <div class="col-6">
            <h6>watchers</h6>
            <p><a href="<?php echo $pipeline->html_url; ?>/watchers"><?php echo $pipeline_metrics[
    'watchers_count'
]; ?></a></p>
        </div>
    </div>

    <div class="row border-bottom">
        <div class="col-6">
            <h6>last release</h6>
            <p><a href="/<?php echo $pipeline->name; ?>/releases"><?php echo $last_release_time; ?></a></p>
        </div>
        <div class="col-6">
            <h6>last updated</h6>
            <p><?php echo $last_commit; ?></p>
        </div>
    </div>

    <div class="row border-bottom">
        <div class="col-6">
            <h6>open issues</h6>
            <p><a href="<?php echo $pipeline->html_url; ?>/issues"><?php echo $pipeline_metrics[
    'open_issues_count'
]; ?></a></p>
        </div>
        <div class="col-6">
            <h6>open pull requests</h6>
            <p><a href="<?php echo $pipeline->html_url; ?>/pulls"><?php echo $pipeline_metrics[
    'open_pr_count'
]; ?></a></p>
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
        <p><a class="btn btn-sm btn-outline-info" href="https://nfcore.slack.com/channels/<?php echo $pipeline->name; ?>"><i class="fab fa-slack me-1"></i> Ask a question on Slack</a></p>
        <p><a class="btn btn-sm btn-outline-secondary" href="<?php echo $pipeline->html_url; ?>/issues"><i class="fab fa-github me-1"></i> Open an issue on GitHub</a></p>
    </div>
</div>

<?php // Collect this content into a variable to be inserted in to the very end of the HTML

ob_start(); ?>

<div class="toast cmd_copied" id="pipeline_sidebar_cmd_copied" role="alert" aria-live="assertive" aria-atomic="true">
    <div class="toast-header">
        <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded me-2" alt="">
        <strong class="me-auto">Command copied!</strong>
        <button type="button" class="ms-2 mb-1 btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
    <div class="toast-body">
        Paste this command into your terminal to run the pipeline with a small test dataset.
    </div>
</div>

<script type="text/javascript">
    $(function() {
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
                datasets: [{
                    backgroundColor: 'rgba(84, 171, 106, 0.2)',
                    borderColor: 'rgba(84, 171, 106, 1)',
                    pointRadius: 0,
                    pointHoverBorderColor: 'rgba(84, 171, 106, 0)', // transparent
                    pointHoverBackgroundColor: 'rgba(84, 171, 106, 0)', // transparent
                    data: [
                        <?php
                        $dates = [];
                        foreach ($traffic_stats as $clone) {
                            if (!is_numeric($clone['clones'])) {
                                continue;
                            }
                            echo '{ x: "' .
                                date('Y-m-d', strtotime($clone['timestamp'])) .
                                '", y: ' .
                                $clone['clones'] .
                                ' },' .
                                "\n\t\t\t";
                        }
                        ?>
                    ]
                }],
            },
            type: 'LineWithLine',
            options: {
                onClick: function(e) {
                    window.location.href = '/<?php echo $pipeline->name; ?>/releases_stats';
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
                            $('#clones_header').text('clones - ' + tooltipItems[0]['label']);
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

