<?php

require_once('../includes/functions.php');
usort($pipeline->releases, 'rsort_releases');

$title = 'nf-core/<br class="d-sm-none">'.$pipeline->name;
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

# Try to fetch the nextflow_schema.json file for the latest release, to see whether we can have a 'Launch' button
$release = 'dev';
if(count($pipeline->releases) > 0){
  $release = $pipeline->releases[0]->tag_name;
}
$gh_launch_schema_fn = dirname(dirname(__FILE__))."/api_cache/json_schema/{$pipeline->name}/{$release}.json";
$gh_launch_no_schema_fn = dirname(dirname(__FILE__))."/api_cache/json_schema/{$pipeline->name}/{$release}.NO_SCHEMA";
# Build directories if needed
if (!is_dir(dirname($gh_launch_schema_fn))) {
  mkdir(dirname($gh_launch_schema_fn), 0777, true);
}
// Load cache if not 'dev'
if((!file_exists($gh_launch_schema_fn) && !file_exists($gh_launch_no_schema_fn)) || $release == 'dev'){
  $api_opts = stream_context_create([ 'http' => [ 'method' => 'GET', 'header' => [ 'User-Agent: PHP' ] ] ]);
  $gh_launch_schema_url = "https://api.github.com/repos/nf-core/{$pipeline->name}/contents/nextflow_schema.json?ref={$release}";
  $gh_launch_schema_json = file_get_contents($gh_launch_schema_url, false, $api_opts);
  if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
    # Remember for next time
    file_put_contents($gh_launch_no_schema_fn, '');
  } else {
    # Save cache
    file_put_contents($gh_launch_schema_fn, $gh_launch_schema_json);
  }
}

# Get details for the Call To Action button
$pipeline_warning = '';
if(count($pipeline->releases) > 0){
  if(file_exists($gh_launch_schema_fn)){
    $cta_btn = '<a href="/launch?pipeline='.$pipeline->name.'&release='.$pipeline->releases[0]->tag_name.'" class="btn btn-success btn-lg"><i class="fad fa-rocket-launch mr-1"></i> Launch version '.$pipeline->releases[0]->tag_name.'</a>';
  } else {
    $cta_btn = '<a href="'.$pipeline->releases[0]->html_url.'" class="btn btn-success btn-lg">Get version '.$pipeline->releases[0]->tag_name.'</a>';
  }
} else {
  if(file_exists($gh_launch_schema_fn)){
    $cta_btn = '<a href="/launch?pipeline='.$pipeline->name.'&release=dev" class="btn btn-success btn-lg"><i class="fad fa-rocket-launch mr-1"></i> Launch development version</a>';
  } else {
    $cta_btn = '<a href="'.$pipeline->html_url.'" class="btn btn-success btn-lg">See the development code</a>';
  }
  $pipeline_warning = '<div class="alert alert-danger">This pipeline is currently in development and does not yet have any stable releases.</div>';
}
if($pipeline->archived){
  $pipeline_warning = '<div class="alert alert-warning">This pipeline has been archived and is no longer being actively maintained.</div>';
}

# Extra HTML for the header - tags and GitHub URL
?>

<div class="mainpage-subheader-heading">
  <div class="container text-center">
    <?php echo $pipeline_warning; ?>
    <p><?php echo $cta_btn; ?></p>
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
    <a class="nav-link<?php if(count($path_parts) > 1 && $path_parts[1] == 'docs'){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>/docs">Doc<span class="d-none d-sm-inline">umentation</span><span class="d-sm-none">s</span></a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($_GET['path'] == $pipeline->name.'/stats'){ echo ' active'; } ?>" href="/<?php echo $pipeline->name; ?>/stats">Stat<span class="d-none d-sm-inline">istic</span>s</a>
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
if($pagetab == 'home' || $pagetab == 'docs' || $pagetab == 'releases'){
  echo '<div class="rendered-markdown">'.$content.'</div>';
} else {
  echo $content;
}

if($pagetab !== 'stats'){
  echo '</div></div>';
}

include('../includes/footer.php');
