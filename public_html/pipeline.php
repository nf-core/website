<?php

require_once('../includes/functions.php');
usort($pipeline->releases, 'rsort_releases');

$title = 'nf-core/<br class="d-sm-none">'.$pipeline->name;
$subtitle = $pipeline->description;

require_once('../includes/pipeline_page/components.php');

# Markdown file - Readme or docs
$pagetab = false;
$release = 'dev';
$release_href = $pipeline->name;
if(count($pipeline->releases) > 0){
  $release = $pipeline->releases[0]->tag_name;
  $release_url = $pipeline->releases[0]->html_url;
}
$schema_content = '';
if(count($path_parts) == 1){
    $pagetab = 'home';
    require_once('../includes/pipeline_page/docs.php');
    $docs = true;
}
# Base readme - redirect to pipeline root
else if($_GET['path'] == $pipeline->name.'/README.md' || $_GET['path'] == $pipeline->name.'/README'){
    header('Location: /'.$pipeline->name);
    exit;
}
# Usage docs
else if($_GET['path'] == $pipeline->name.'/usage'){
    $pagetab = 'usage';
    require_once('../includes/pipeline_page/usage.php');
}
# output docs
else if($_GET['path'] == $pipeline->name.'/output'){
    $pagetab = 'output';
    require_once('../includes/pipeline_page/output.php');
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
// handling urls with a version in it
else if(count($path_parts) >= 2 && preg_match('/^\d+\.\d+\.\d+|^\d+\.\d+$|^dev/',$path_parts[1])){
  $release = $path_parts[1];
  foreach($pipeline->releases as $releases)
  if($releases==$release){
    $release_url = $releases->html_url;
  }
  $release_href = $pipeline->name.'/'.$release;
  if(count($path_parts) == 2 && preg_match('/^\d+\.\d+\.\d+|^\d+.\d+$|^dev/',$path_parts[1])){
    $pagetab = 'home';
    require_once('../includes/pipeline_page/docs.php');
  }else{

  switch($path_parts[2]){
    case 'usage':
      $pagetab = 'usage';
      require_once('../includes/pipeline_page/usage.php');
    break;
    case 'output':
      $pagetab = 'output';
      require_once('../includes/pipeline_page/output.php');
      break;
    case 'stats':
      $pagetab = 'stats';
      require_once('../includes/pipeline_page/stats.php');
      break;
    case 'releases':
      $pagetab = 'releases';
      require_once('../includes/pipeline_page/releases.php');
      break;
    default:
      header('HTTP/1.1 404 Not Found');
      include('404.php');
      die();
  }


  }

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
  if($release!="dev"){
  if(file_exists($gh_launch_schema_fn)){
    $cta_btn = '<a href="/launch?pipeline='.$pipeline->name.'&release='.$release.'" class="btn btn-success btn-lg"><i class="fad fa-rocket-launch mr-1"></i> Launch version '.$release.'</a>';
  } else {
    $cta_btn = '<a href="'.$release_url.'" class="btn btn-success btn-lg"><i class="fas fa-tags mr-1"></i> See version '.$release.'</a>';
  }
  }else{
    if(file_exists($gh_launch_schema_fn)){
    $cta_btn = '<a href="/launch?pipeline='.$pipeline->name.'&release=dev" class="btn btn-success btn-lg"><i class="fad fa-rocket-launch mr-1"></i> Launch development version</a>';
  } else {
    $cta_btn = '<a href="'.$pipeline->html_url.'" class="btn btn-success btn-lg"><i class="fad fa-construction mr-1"></i> See the latest code</a>';
  }
  }
} else {
  if(file_exists($gh_launch_schema_fn)){
    $cta_btn = '<a href="/launch?pipeline='.$pipeline->name.'&release=dev" class="btn btn-success btn-lg"><i class="fad fa-rocket-launch mr-1"></i> Launch development version</a>';
  } else {
    $cta_btn = '<a href="'.$pipeline->html_url.'" class="btn btn-success btn-lg"><i class="fad fa-construction mr-1"></i> See the latest code</a>';
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
    <a class="nav-link<?php if(count($path_parts) == 1){ echo ' active'; } ?>" href="/<?php echo $release_href; ?>">Readme</a>
  </li>

  <li class="nav-item">
    <a class="nav-link<?php if($pagetab=='usage'){ echo ' active'; } ?>" href="/<?php echo $release_href; ?>/usage">Usage</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($pagetab=='output'){ echo ' active'; } ?>" href="/<?php echo $release_href; ?>/output">Outputs</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($pagetab=='stats'){ echo ' active'; } ?>" href="/<?php echo $release_href; ?>/stats">Stat<span class="d-none d-sm-inline">istic</span>s</a>
  </li>
  <li class="nav-item">
    <a class="nav-link<?php if($pagetab=='releases'){ echo ' active'; } ?>" href="/<?php echo $release_href; ?>/releases">Releases</a>
  </li>
  <?php if($pagetab == 'home' || $pagetab == 'output' || $pagetab == 'usage'): ?>
  <li>
    <div class="input-group">
      <div class="input-group-prepend">
        <label class="input-group-text" for="version_select"><i class="fas fa-tags"></i></label>
      </div>
      <select class="custom-select" id="version_select" data-pipeline="<?php echo $pipeline->name?>">
        <?php
          foreach($pipeline->releases as $release){
          ?>
            <option value="<?php echo strtolower($release->tag_name); ?>"><?php echo $release->tag_name; ?></option>
          <?php
          }
          ?>
          <option value="dev">dev</option>
      </select>
    </div>
  </li>
  <?php endif; ?>
</ul>

<?php

# echo '<pre>'.print_r($pipeline, true).'</pre>';
# echo '<pre>'.print_r($gh_tree_json, true).'</pre>';

if($pagetab !== 'stats'){
    echo '<div class="row"><div class="col-lg-8 order-lg-1">';
}
# Print content
if($pagetab == 'home' || $pagetab == 'output' || $pagetab == 'usage' || $pagetab == 'releases'){
  if(preg_match('/<!-- params-docs -->/')){
    $content = '<div class="rendered-markdown">'.preg_replace('/<!-- params-docs -->/',$schema_content,$content).'</div>';
  } else {
    $content = '<div class="rendered-markdown">'.$content.$schema_content.'</div>';
  }
  echo $content;
}
else {
  echo $content;
}
if($pagetab !== 'stats'){
    echo '</div><div class="col-lg-4 order-lg-12"><div class="side-sub-subnav sticky-top">';
    if($pagetab == 'usage'){
      $toc = '<div class="btn-group mt-1" role="group">
                    <button class="btn btn-outline-secondary collapse-groups-btn" id="toggle_details" data-toggle="collapse"  data-target=".schema-docs-help-text" aria-expanded="false"><i class="fa mr-1"></i> Toggle details</button>
                    <button class="btn btn-outline-secondary collapse-groups-btn" id="show_hidden" data-toggle="collapse" data-target=".param_hidden" aria-expanded="false"><i class="fa mr-1"></i> Show hidden</button>
                </div>';
      $toc .= '<nav class="toc nav flex-column">'.generate_toc($content).'</nav>';
        echo $toc;
    } else if($pagetab == 'output'){
      $toc = '<nav class="toc nav flex-column">'.generate_toc($content).'</nav>';
        echo $toc;
    } else {
        echo $pipeline_stats_sidebar;
    }
}

if($pagetab !== 'stats'){
  echo '</div></div></div>';
}

include('../includes/footer.php');
