<?php

require_once('../includes/functions.php');
usort($pipeline->releases, 'rsort_releases');

$title = '<a href="/'.$pipeline->name.'">nf-core/<br class="d-sm-none">'.$pipeline->name.'</a>';
$subtitle = $pipeline->description;

require_once('../includes/pipeline_page/components.php');

# Markdown file - Readme or docs
$pagetab = false;
if(count($path_parts) == 1){
    $pagetab = 'home';
    require_once('../includes/pipeline_page/docs.php');
    $docs = true;
} else if($path_parts[1] == 'docs'){
    $pagetab = 'docs';
    require_once('../includes/pipeline_page/docs.php');
    $docs = true;
}
# Base readme - redirect to pipeline root
else if($_GET['path'] == $pipeline->name.'/README.md' || $_GET['path'] == $pipeline->name.'/README'){
    header('Location: /'.$pipeline->name);
    exit;
}
# Stats
else if($_GET['path'] == $pipeline->name.'/stats'){
    $pagetab = 'stats';
    require_once('../includes/pipeline_page/stats.php');
}
# Releases
else if($_GET['path'] == $pipeline->name.'/releases'){
    $pagetab = 'releases';
    require_once('../includes/pipeline_page/releases.php');
}
# Some other URL pattern that we don't recognise - 404
else {
    header('HTTP/1.1 404 Not Found');
    include('404.php');
    die();
}


# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
include('../includes/header.php');

# Pipeline subheader

# Get details for the button to the latest release
$dev_warning = '';
$archived_warning = '';
if(count($pipeline->releases) > 0){
    $download_bn = '<a href="'.$pipeline->releases[0]->html_url.'" class="btn btn-success btn-lg">Get version '.$pipeline->releases[0]->tag_name.'</a>';
} else {
    $download_bn = '<a href="'.$pipeline->html_url.'" class="btn btn-success btn-lg">See the development code</a>';
    $dev_warning = '<div class="alert alert-danger">This pipeline is currently in development and does not yet have any stable releases.</div>';
}
if($pipeline->archived){
  $archived_warning = '<div class="alert alert-warning">This pipeline has been archived and is no longer being actively maintained.</div>';
}

# Extra HTML for the header - tags and GitHub URL
?>

<div class="mainpage-subheader-heading">
  <div class="container text-center">
    <?php echo $dev_warning.$archived_warning; ?>
    <p><?php echo $download_bn; ?></p>
    <p class="mb-0"><a href="<?php echo $pipeline->html_url; ?>" class="text-dark"><i class="fab fa-github"></i> <?php echo $pipeline->html_url; ?></a></p>
  </div>
</div>
<div class="triangle subheader-triangle-down"></div>

<div class="container main-content">

<ul class="nav nav-fill nfcore-subnav">
  <li class="nav-item">
    <a class="nav-link<?php if(count($path_parts) == 1){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>">Readme</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if(count($path_parts) > 1 && $path_parts[1] == 'docs'){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>/docs">Documentation</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($_GET['path'] == $pipeline->name.'/stats'){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>/stats">Stats</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($_GET['path'] == $pipeline->name.'/releases'){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>/releases">Releases</a>
  </li>
</ul>

<?php

# echo '<pre>'.print_r($pipeline, true).'</pre>';
# echo '<pre>'.print_r($gh_tree_json, true).'</pre>';

if($pagetab !== 'stats'){
    echo '<div class="row"><div class="col-lg-4 order-lg-12"><div class="side-sub-subnav sticky-top">';
    if($pagetab == 'docs'){
        echo '<div class="pipeline-page-toc">'.$md_toc_html.'</div>';
    } else {
        echo $pipeline_stats_sidebar;
    }
    echo '</div></div><div class="col-lg-8 order-lg-1">';
}

# Print content
if($pagetab == 'home' || $pagetab == 'docs'){
  echo '<div class="rendered-markdown">'.$content.'</div>';
} else {
  echo $content;
}

if($pagetab !== 'stats'){
  echo '</div></div>';
}

include('../includes/footer.php');
