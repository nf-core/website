<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

// Grab the libraries for the markdown parsing
require_once(dirname(dirname(__FILE__)).'/libraries/parsedown/Parsedown.php');
require_once(dirname(dirname(__FILE__)).'/libraries/parsedown-extra/ParsedownExtra.php');

ob_start();
foreach($pipeline->releases as $release){

    $pd = new ParsedownExtra();
    $body_html = $pd->text($release->body);

    ?>

<h1 id="<?php echo $release->tag_name; ?>"><a href="#<?php echo $release->tag_name; ?>" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
  <a href="<?php echo $release->html_url; ?>" style="color: inherit;">
    <span class="badge badge-secondary float-right"><i class="fas fa-tag"></i> <?php echo $release->tag_name; ?></span>
    <?php echo $release->name; ?>
  </a>
</h1>
<p>
  <span class="float-right">
    <a href="<?php echo $release->zipball_url; ?>" class="btn btn-sm btn-outline-success">Download .zip</a>
    <a href="<?php echo $release->tarball_url; ?>" class="btn btn-sm btn-outline-success">Download .tar.gz</a>
  </span>
  <i class="far fa-calendar-alt"></i> Released <?php echo date('j M Y', strtotime($release->published_at)); ?>
  <small class="text-muted">(<?php echo time_ago($release->published_at); ?>)</small>
</p>

<p><?php echo $body_html; ?></p>

<?php
}

$content = ob_get_contents();
ob_end_clean();
