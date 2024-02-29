<?php

$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;

$msg = '';
if (isset($_GET['update']) && time() - $pipelines_json->updated > 60 * 10) {
    $output = shell_exec('php ../update_pipeline_details.php 2>&1 | tee -a /home/nfcore/update.log');
    $msg = '<div class="alert alert-success">Manual pipeline update sync triggered: <pre>' . $output . '</pre></div>';
} elseif (isset($_GET['update'])) {
    $msg =
        '<div class="alert alert-warning">Manual pipeline update can only be triggered at least 10 minutes apart</div>';
}

$title = 'Pipelines';
$subtitle =
    'Browse the <strong>' .
    $pipelines_json->pipeline_count .
    '</strong> pipelines that are currently available as part of sanger-tol.';
include '../includes/header.php';

usort($pipelines, 'rsort_pipelines');

echo $msg;
?>

<h1>Available Pipelines</h1>
<p class="mb-3">Can you think of another pipeline that would fit in well?
  <a href="https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato%2Cmuffato&labels=connect&projects=&template=contact_us.yaml&title=%5BContact+Us%5D%3A+">Let us know!</a>
</p>

<div class="btn-toolbar mb-4 pipelines-toolbar" role="toolbar">
  <div class="pipeline-filters input-group input-group-sm me-2 mt-2">
    <input type="search" class="form-control" placeholder="Search keywords" value="<?php echo isset($_GET['q'])
        ? htmlspecialchars($_GET['q'], ENT_QUOTES, 'UTF-8')
        : ''; ?>">
  </div>
  <div class="btn-group btn-group-sm mt-2 d-none d-lg-block" role="group">
    <button type="button" class="btn text-body">Filter:</button>
  </div>
  <div class="pipeline-filters btn-group btn-group-sm me-2 mt-2">
    <?php
    if ($pipelines_json->published_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-bs-target=".pipeline-released">
        Released
        <span class="badge bg-light text-dark d-none d-lg-inline"><?php echo $pipelines_json->published_count; ?></span>
      </button>
    <?php endif;
    if ($pipelines_json->devel_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-bs-target=".pipeline-dev">
        <span class="d-none d-lg-inline">Under development</span> <span class="d-inline d-lg-none">Dev</span>
        <span class="badge bg-light text-dark d-none d-lg-inline"><?php echo $pipelines_json->devel_count; ?></span>
      </button>
    <?php endif;
    if ($pipelines_json->archived_count > 0): ?>
      <button type="button" class="btn btn-sm btn-outline-success active" data-bs-target=".pipeline-archived">
        Archived
        <span class="badge bg-light text-dark d-none d-lg-inline"><?php echo $pipelines_json->archived_count; ?></span>
      </button>
    <?php endif;
    ?>
  </div>
  <div class="btn-group btn-group-sm mt-2 d-none d-lg-block" role="group">
    <button type="button" class="btn text-body">Sort:</button>
  </div>
  <div class="pipeline-sorts btn-group btn-group-sm me-2 mt-2" role="group">
    <button type="button" class="btn btn-outline-success active"><span class="d-none d-xl-inline">Last</span> Release</button>
    <button type="button" class="btn btn-outline-success">Alphabetical</button>
    <button type="button" class="btn btn-outline-success">Stars</button>
  </div>
  <div class="btn-group btn-group-sm mt-2 d-none d-xl-block" role="group">
    <button type="button" class="btn text-body">Display:</button>
  </div>
  <div class="btn-group btn-group-sm mt-2" role="group">
    <button data-dtype="blocks" type="button" class="display-btn btn btn-outline-success active" title="Display as blocks" data-bs-toggle="tooltip"><i class="fas fa-th-large"></i></button>
    <button data-dtype="list" type="button" class="display-btn btn btn-outline-success" title="Display as list" data-bs-toggle="tooltip"><i class="fas fa-bars"></i></button>
  </div>
</div>

<p class="no-pipelines text-muted mt-5" style="display: none;">No pipelines found..</p>

<div class="row row-cols-1 row-cols-md-2 g-4 pipelines-container">
  <?php foreach ($pipelines as $wf): ?>
    <div class="col">
      <div class="card h-100 pipeline <?php if ($wf->archived): ?>pipeline-archived<?php elseif (
          count($wf->releases) == 0
      ): ?>pipeline-dev<?php else: ?>pipeline-released<?php endif; ?>">
        <div class="card-body clearfix">
          <h3 class="card-title mb-0">
            <?php if ($wf->stargazers_count > 0): ?>
              <a href="<?php echo $wf->html_url; ?>/stargazers" target="_blank" class="stargazers text-decoration-none mt-2 ms-2" title="<?php echo $wf->stargazers_count; ?> stargazers on GitHub <small class='fas fa-external-link-alt ms-2'></small>" data-bs-toggle="tooltip" data-html="true">
                <i class="far fa-star"></i>
                <?php echo $wf->stargazers_count; ?>
              </a>
            <?php endif; ?>
            <a href="/<?php echo $wf->name; ?>" class="pipeline-name">
              <span class="d-none d-lg-inline">sanger-tol/</span><?php echo $wf->name; ?>
            </a>
            <?php if ($wf->archived): ?>
              <small class="status-icon text-warning ms-2 fas fa-archive" title="This pipeline has been archived and is no longer being maintained." data-bs-toggle="tooltip"></small>
            <?php elseif (count($wf->releases) == 0): ?>
              <small class="status-icon text-danger ms-2 fas fa-wrench" title="This pipeline is under active development. Once released on GitHub, it will be production-ready." data-bs-toggle="tooltip"></small>
            <?php else: ?>
              <small class="status-icon text-success ms-2 fas fa-check" title="This pipeline is released, tested and good to go." data-bs-toggle="tooltip"></small>
            <?php endif; ?>
          </h3>
          <?php if (count($wf->topics) > 0): ?>
            <p class="topics mb-0">
              <?php foreach ($wf->topics as $topic): ?>
                <a href="/pipelines?q=<?php echo $topic; ?>" class="badge pipeline-topic"><?php echo $topic; ?></a>
              <?php endforeach; ?>
            </p>
          <?php endif; ?>
          <p class="card-text mb-0 mt-2"><?php echo $wf->description; ?></p>
          <p class="mb-0 mt-2 dl-btn-row">
            <?php if (count($wf->releases) > 0):
                usort($wf->releases, 'rsort_releases'); ?>
              <a href="<?php echo $wf->releases[0]
                  ->html_url; ?>" target="_blank" class="btn btn-sm btn-outline-success">
                Version <strong><?php echo $wf->releases[0]->tag_name; ?></strong>
              </a> &nbsp;
              <small class="text-black-50 publish-date" data-pubdate="<?php echo strtotime(
                  $wf->releases[0]->published_at,
              ); ?>">Published <?php echo time_ago($wf->releases[0]->published_at); ?></small>
            <?php
            else:
                 ?>
              <small class="text-danger">No releases yet</small>
            <?php
            endif; ?>
          </p>
        </div>
      </div>
    </div>
  <?php endforeach; ?>
</div>

<p class="mt-5 small text-muted">
  Page last synced with GitHub <?php echo time_ago($pipelines_json->updated); ?>.
  <a href="/pipelines?update">Click here</a> to trigger an update</a>.
  See also <a href="/pipeline_health">pipeline repository health</a>.
</p>

<?php include '../includes/footer.php'; ?>
