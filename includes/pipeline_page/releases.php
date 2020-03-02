<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

ob_start();

echo '<h1>Version history</h1>';

foreach($pipeline->releases as $release){ ?>

<div class="row">
  <div class="col-auto">
    <a href="#download-<?php echo $release->tag_sha; ?>" class="text-body" data-toggle="collapse">
      <samp><?php echo $release->tag_name; ?></samp>
    </a>
  </div>
  <div class="col">
  </div>
  <div class="col-auto">
    <a href="#download-<?php echo $release->tag_sha; ?>" class="text-body" data-toggle="collapse"><small class="text-muted"><?php echo time_ago($release->published_at); ?></small></a>
    <button class="btn btn-sm btn-link text-body" type="button" data-toggle="collapse" data-target="#download-<?php echo $release->tag_sha; ?>">
      <i class="fas fa-caret-left"></i>
  </button>
  </div>
</div>
<div class="collapse" id="download-<?php echo $release->tag_sha; ?>">
  <div class="row pb-2">
    <div class="col-sm-6 small">
      Released <?php echo date('j M Y', strtotime($release->published_at)); ?> &mdash;
      <code><?php echo substr($release->tag_sha, 0, 7); ?></code>
    </div>
    <div class="col-sm-6 text-right">
      <a href="<?php echo $release->zipball_url; ?>" class="btn btn-sm btn-outline-success">Download .zip</a>
      <a href="<?php echo $release->tarball_url; ?>" class="btn btn-sm btn-outline-success">Download .tar.gz</a>
      <a href="<?php echo $release->html_url; ?>" class="btn btn-sm btn-success">View release</a>
    </div>
  </div>
</div>
<hr class="m-0">



<?php
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

<?php

$content = ob_get_contents();
ob_end_clean();
