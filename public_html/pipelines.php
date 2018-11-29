<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

function rsort_releases($a, $b){
    $t1 = strtotime($a->published_at);
    $t2 = strtotime($b->published_at);
    return $t2 - $t1;
}
function rsort_pipelines($a, $b){
    $t1 = strtotime($a->last_release);
    $t2 = strtotime($b->last_release);
    return $t2 - $t1;
}

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
usort($pipelines, 'rsort_pipelines');

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
$subtitle = 'Browse the <strong>'.$pipelines_json->pipeline_count.'</strong> pipelines that are currently available as part of nf-core.';
include('../includes/header.php');
?>

<h1>Available Pipelines</h1>
<p class="mb-3">Can you think of another pipeline that would fit in well?
<a href="https://nf-core-invite.herokuapp.com/">Let us know!</a></p>

<div class="btn-toolbar mb-4 pipelines-toolbar" role="toolbar">
  <div class="pipeline-filters input-group input-group-sm mr-2 mt-2">
    <input type="text" class="form-control" placeholder="Search keywords">
  </div>
  <div class="btn-group btn-group-sm mt-2 d-none d-sm-block" role="group">
    <button type="button" class="btn btn-link text-body">Filter:</button>
  </div>
  <div class="pipeline-filters btn-group btn-group-sm mr-2 mt-2">
    <?php if($pipelines_json->published_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-target=".pipeline-released">Released <span class="badge badge-light"><?php echo $pipelines_json->published_count; ?></span></button>
    <?php endif;
    if($pipelines_json->devel_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-target=".pipeline-dev">Under development <span class="badge badge-light"><?php echo $pipelines_json->devel_count; ?></span></button>
    <?php endif;
    if($pipelines_json->archived_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-target=".pipeline-archived">Archived <span class="badge badge-light"><?php echo $pipelines_json->archived_count; ?></span></button>
    <?php endif; ?>
  </div>
  <div class="btn-group btn-group-sm mt-2 d-none d-sm-block" role="group">
    <button type="button" class="btn btn-link text-body">Sort:</button>
  </div>
  <div class="pipeline-sorts btn-group btn-group-sm mr-2 mt-2" role="group">
    <button type="button" class="btn btn-outline-success active">Last Release</button>
    <button type="button" class="btn btn-outline-success">Alphabetical</button>
    <button type="button" class="btn btn-outline-success">Status</button>
    <button type="button" class="btn btn-outline-success">Stars</button>
  </div>
</div>

<p class="no-pipelines text-muted mt-5" style="display: none;">No pipelines found..</p>

<div class="card-deck pipelines-container">
<?php foreach($pipelines as $wf): ?>
    <div class="card card_deck_card pipeline <?php if($wf->archived): ?>pipeline-archived<?php elseif(count($wf->releases) == 0): ?>pipeline-dev<?php else: ?>pipeline-released<?php endif; ?>">
        <div class="card-body">
            <h3 class="card-title mb-0">
                <?php if($wf->stargazers_count > 0): ?>
                <a href="<?php echo $wf->html_url; ?>/stargazers" target="_blank" class="stargazers mt-2 ml-2" title="<?php echo $wf->stargazers_count; ?> stargazers on GitHub <small class='fas fa-external-link-alt ml-2'></small>" data-toggle="tooltip" data-html="true">
                    <i class="far fa-star"></i>
                    <?php echo $wf->stargazers_count; ?>
                </a>
                <?php endif; ?>
                <a href="<?php echo $wf->html_url; ?>" target="_blank" class="pipeline-name">
                    <?php echo $wf->full_name; ?>
                </a>
                <?php if($wf->archived): ?>
                <small class="status-icon text-warning ml-2 fas fa-archive" title="This pipeline has been archived and is no longer being maintained." data-toggle="tooltip"></small>
                <?php elseif(count($wf->releases) == 0): ?>
                    <small class="status-icon text-danger ml-2 fas fa-exclamation-triangle" title="This pipeline is under active development. Once released on GitHub, it will be production-ready." data-toggle="tooltip"></small>
                <?php else: ?>
                    <small class="status-icon text-success ml-2 fas fa-check" title="This pipeline is released, tested and good to go." data-toggle="tooltip"></small>
                <?php endif; ?>
            </h3>
            <?php if(count($wf->topics) > 0): ?>
              <p class="topics mb-0">
              <?php foreach($wf->topics as $topic): ?>
                <span class="badge pipeline-topic"><?php echo $topic; ?></span>
              <?php endforeach; ?>
              </p>
            <?php endif; ?>
            <p class="card-text mb-0 mt-2"><?php echo $wf->description; ?></p>
            <p class="mb-0 mt-2">
            <?php if(count($wf->releases) > 0):
                usort($wf->releases, 'rsort_releases');
                ?>
                <a href="<?php echo $wf->releases[0]->html_url; ?>" target="_blank"  class="btn btn-sm btn-outline-success">
                    Version <strong><?php echo $wf->releases[0]->tag_name; ?></strong>
                </a> &nbsp;
                <small class="text-black-50 publish-date" data-pubdate="<?php echo strtotime($wf->releases[0]->published_at); ?>">Published <?php echo time_ago($wf->releases[0]->published_at); ?></small>
            <?php else: ?>
                <small class="text-danger">No releases yet</small>
            <?php endif; ?>
            </p>
        </div>
    </div>
<?php endforeach; ?>
</div>

<p class="mt-5"><small class="text-muted">Page last synced with GitHub <?php echo time_ago($pipelines_json->updated); ?>.</small></p>

<?php include('../includes/footer.php'); ?>
