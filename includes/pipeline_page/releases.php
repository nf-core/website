<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

ob_start();

echo '<h1>Version history</h1>';

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
      <a href="<?php echo $releases->html_url; ?>" class="btn btn-sm btn-success">View release</a>
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

<?php

$content = ob_get_contents();
ob_end_clean();
