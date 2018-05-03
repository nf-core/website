<?php

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

// From https://stackoverflow.com/a/18891474/713980
function time_ago($date) {
    $periods = array("second", "minute", "hour", "day", "week", "month", "year", "decade");
    $lengths = array("60", "60", "24", "7", "4.35", "12", "10");
    $now = time();
    if(is_numeric($date)) $unix_date = $date;
    else $unix_date = strtotime($date);
    // check validity of date
    if (empty($unix_date)) {
        return $date;
    }
    // is it future date or past date
    if ($now > $unix_date) {
        $difference = $now - $unix_date;
        $tense = "ago";
    } else {
        $difference = $unix_date - $now;
        $tense = "from now";
    }
    for ($j = 0; $difference >= $lengths[$j] && $j < count($lengths) - 1; $j++) {
        $difference /= $lengths[$j];
    }
    $difference = round($difference);
    if ($difference != 1) {
        $periods[$j].= "s";
    }
    return "$difference $periods[$j] {$tense}";
}

$title = 'Pipelines';
$subtitle = 'Browse the '.$pipelines_json->pipeline_count.' pipelines that are currently available as part of nf-core.';
include('../includes/header.php');
?>

<h1>Available Pipelines</h1>
<p>Can you think of another pipeline that would fit in well?
<a href="https://github.com/nf-core/nf-core.github.io/issues/new">Let us know!</a></p>

<?php if($pipelines_json->published_count > 0): ?>

<h2>Production Pipelines</h2>

<div class="row">
<?php foreach($pipelines as $wf): if(count($wf->releases) > 0 and !$wf->archived): ?>
    <div class="col-sm-6 mb-3">
        <div class="card pipeline">
            <div class="card-body">
                <h3 class="card-title">
                    <a href="<?php echo $wf->html_url; ?>" target="_blank" class="stargazers mt-2">
                        <i class="far fa-star"></i>
                        <?php echo $wf->stargazers_count; ?>
                    <a href="<?php echo $wf->html_url; ?>" target="_blank">
                        <?php echo $wf->full_name; ?>
                    </a>
                </h3>
                <p class="card-text"><?php echo $wf->description; ?></p>
                <p class="mb-0">
                    <a href="<?php echo $wf->releases[0]->html_url; ?>" target="_blank"  class="btn btn-sm btn-outline-success">
                        Version <strong><?php echo $wf->releases[0]->tag_name; ?></strong>
                    </a> &nbsp;
                    <small class="text-black-50">Published <?php echo time_ago($wf->releases[0]->published_at); ?></small>
                </p>
            </div>
        </div>
    </div>
<?php endif; endforeach; ?>
</div>

<?php endif;

if($pipelines_json->devel_count > 0): ?>

<h2>Pipelines in development</h2>
<p>These pipelines aren't yet ready for the prime-time, but are under active development.
Once they've had a release on GitHub, they will be classed as production-ready pipelines.</p>

<div class="row">
<?php foreach($pipelines as $wf): if(count($wf->releases) == 0): ?>
    <div class="col-sm-6 mb-3">
        <div class="card pipeline">
            <div class="card-body">
                <h3 class="card-title">
                    <a href="<?php echo $wf->html_url; ?>" target="_blank" class="stargazers mt-2">
                        <i class="far fa-star"></i>
                        <?php echo $wf->stargazers_count; ?>
                    <a href="<?php echo $wf->html_url; ?>" target="_blank">
                        <?php echo $wf->full_name; ?>
                    </a>
                </h3>
                <p class="card-text"><?php echo $wf->description; ?></p>
                <p class="mb-0">
                    <small class="text-danger">No releases yet</small>
                </p>
            </div>
        </div>
    </div>
<?php endif; endforeach; ?>
</div>

<?php endif;

if($pipelines_json->archived_count > 0): ?>

<h2>Archived Pipelines</h2>
<p>These pipelines have been archived and are no longer under active development.</p>

<div class="row">
<?php foreach($pipelines as $wf): if($wf->archived): ?>
    <div class="col-sm-6 mb-3">
        <div class="card pipeline">
            <div class="card-body">
                <h3 class="card-title">
                    <a href="<?php echo $wf->html_url; ?>" target="_blank" class="stargazers mt-2">
                        <i class="far fa-star"></i>
                        <?php echo $wf->stargazers_count; ?>
                    <a href="<?php echo $wf->html_url; ?>" target="_blank">
                        <?php echo $wf->full_name; ?>
                    </a>
                </h3>
                <p class="card-text"><?php echo $wf->description; ?></p>
                <p class="mb-0">
                    <a href="<?php echo $wf->releases[0]->html_url; ?>" target="_blank"  class="btn btn-sm btn-outline-success">
                        Version <strong><?php echo $wf->releases[0]->tag_name; ?></strong>
                    </a> &nbsp;
                    <small class="text-black-50">Published <?php echo time_ago($wf->releases[0]->published_at); ?></small>
                </p>
            </div>
        </div>
    </div>
<?php endif; endforeach; ?>
</div>

<?php endif; ?>

<p class="mt-5"><small class="text-muted">Page last synced with GitHub <?php echo time_ago($pipelines_json->updated); ?>.</small></p>

<?php include('../includes/footer.php'); ?>
