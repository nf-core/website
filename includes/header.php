<?php

require_once('functions.php');

// Find the latest commit hash to prevent caching assets
$git_sha = trim(shell_exec("cd ".dirname(__FILE__)." && git rev-parse --short=7 HEAD"));
if(strlen($git_sha) != 7){
  $git_sha = '';
}

// Theme switcher cookie
$theme = 'auto';
if(isset($_COOKIE['nfcoretheme']) && in_array($_COOKIE['nfcoretheme'], ['auto', 'light', 'dark'])) {
  $theme = $_COOKIE['nfcoretheme'];
}

// Convert Markdown to HTML if a filename is given
if(isset($markdown_fn) and $markdown_fn){
  require_once('parse_md.php');
  $parsed_out = parse_md($markdown_fn);
  $content = $parsed_out["content"];
  $meta = $parsed_out["meta"];
  $title = $parsed_out["title"];
  $subtitle = $parsed_out["subtitle"];
}

// Page title
$page_title = 'nf-core';
if(isset($title) && strlen($title) > 0){
  $page_title = preg_replace('/^nf-core\//', '', strip_tags($title)).' &raquo; nf-core';
}
// Page meta description
$page_meta = 'A collection of high quality Nextflow pipelines';
if(isset($subtitle) && strlen($subtitle) > 0){
  $page_meta = strip_tags($subtitle);
}

?><!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title><?php echo $page_title; ?></title>
    <meta name="description" content="<?php echo $page_meta; ?>">
    <meta name="author" content="Phil Ewels">
    <link rel="shortcut icon" href="/assets/img/logo/nf-core-logo-square.png" type="image/png" />
    <link href="/assets/css/bootstrap.min.css" rel="stylesheet">
    <link href="/assets/css/code_highlighting/github.css" rel="stylesheet" >
    <link href="/assets/css/Chart.min.css" rel="stylesheet">
    <!-- FontAwesome -->
    <script src="https://kit.fontawesome.com/471b59d3f8.js"></script>
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-68098153-2"></script>
    <!-- Other JS -->
    <script src="/assets/js/jquery-3.4.1.min.js"></script>
    <script src="/assets/js/popper.min.js"></script>
    <script src="/assets/js/bootstrap.min.js"></script>
    <script src="/assets/js/highlight.pack.js"></script>
    <!-- Page-specific CSS and JS -->
    <?php if(isset($import_leaflet) && $import_leaflet): ?>
    <link href="/assets/css/leaflet.css" rel="stylesheet">
    <link href="/assets/css/leaflet.fullscreen.css" rel="stylesheet">
    <script src="/assets/js/leaflet.js"></script>
    <script src="/assets/js/Leaflet.fullscreen.min.js"></script>
    <?php endif;
    if(isset($import_chartjs) && $import_chartjs): ?>
    <script src="/assets/js/moment.js"></script>
    <script src="/assets/js/Chart.min.js"></script>
    <script src="/assets/js/hammer.min.js"></script>
    <script src="/assets/js/chartjs-plugin-zoom.min.js"></script>
    <script src="/assets/js/canvas2svg.js"></script>
    <script src="/assets/js/FileSaver.js"></script>
    <?php endif;
    if(isset($import_schema_launcher) && $import_schema_launcher): ?>
    <script src="/assets/js/moment.js"></script>
    <script src="/assets/js/showdown.min.js"></script>
    <script src="/assets/js/nf-core-schema-launcher.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif;
    if(isset($aws) && $aws): ?>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.12.0/css/all.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.22/css/dataTables.bootstrap4.min.css">

    <script src="https://cdnjs.cloudflare.com/ajax/libs/bootbox.js/5.4.0/bootbox.min.js"></script>
    <script src="https://sdk.amazonaws.com/js/aws-sdk-2.765.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.22/js/dataTables.bootstrap4.min.js"></script>
    <script src="/assets/js/aws-s3-explorer.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif;
    if(isset($import_schema_builder) && $import_schema_builder): ?>
    <link href="/assets/css/jquery-ui.min.css" rel="stylesheet">
    <script src="/assets/js/jquery-ui.min.js"></script>
    <script src="/assets/js/moment.js"></script>
    <script src="/assets/js/showdown.min.js"></script>
    <script src="/assets/js/nf-core-schema-builder.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif; ?>
    <script src="/assets/js/jquery.tablesorter.min.js"></script>

    <!-- Custom nf-core CSS and JS -->
    <link href="/assets/css/nf-core.css?c=<?php echo $git_sha; ?>" rel="stylesheet">
    <link href="/assets/css/nf-core-<?php echo $theme; ?>.css?c=<?php echo $git_sha; ?>" rel="stylesheet" id="theme-stylesheet">
    <script src="/assets/js/nf-core.js?c=<?php echo $git_sha; ?>"></script>
    <script>window.dataLayer = window.dataLayer || []; function gtag(){dataLayer.push(arguments);}  gtag('js', new Date()); gtag('config', 'UA-68098153-2'); </script>
  </head>
  <body data-spy="scroll" data-target=".toc" data-offset="15">
    <nav class="navbar fixed-top navbar-expand-md navbar-light site-nav">
      <a class="navbar-brand d-md-none" href="/">
        <img height="25px" src="/assets/img/logo/nf-core-logo.svg" class="hide-dark">
        <img height="25px" src="/assets/img/logo/nf-core-logo-darkbg.svg" class="hide-light">
      </a>
      <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarCollapse" aria-controls="navbarCollapse" aria-expanded="false" aria-label="Toggle navigation">
        <span class="navbar-toggler-icon"></span>
      </button>
      <div class="collapse navbar-collapse justify-content-md-center" id="navbarCollapse">
        <ul class="navbar-nav">
          <li class="nav-item p-1">
            <a class="nav-link" href="/">Home</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/pipelines">Pipelines</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/tools">Tools</a>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link" href="/usage/introduction" role="button" data-toggle="dropdown">Usage</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/usage/introduction">Getting started</a>
              <a class="dropdown-item" href="/usage/installation">Installation</a>
              <a class="dropdown-item" href="/usage/configuration">Pipeline configuration</a>
              <a class="dropdown-item" href="/usage/offline">Running offline</a>
              <a class="dropdown-item" href="/usage/nf_core_tutorial">nf-core tutorial</a>
              <a class="dropdown-item" href="/usage/reference_genomes">Reference genomes</a>
              <a class="dropdown-item" href="/usage/data_management">Data Management</a>
              <a class="dropdown-item" href="/usage/troubleshooting">Troubleshooting</a>
              <a class="dropdown-item" href="/usage/nextflow">Nextflow resources</a>
            </div>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link" href="/developers/guidelines" role="button" data-toggle="dropdown">Developers</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/developers/guidelines">Guidelines</a>
              <a class="dropdown-item" href="/developers/adding_pipelines">Adding a new pipeline</a>
              <a class="dropdown-item" href="/developers/release_checklist">Release checklist</a>
              <a class="dropdown-item" href="/errors">Lint error codes</a>
              <a class="dropdown-item" href="/developers/sync">Template synchronisation</a>
              <a class="dropdown-item" href="/developers/design_guidelines">Design Guidelines</a>
            </div>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link" href="/about" role="button" data-toggle="dropdown">About</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/about">About nf-core</a>
              <a class="dropdown-item" href="/events">Events</a>
              <a class="dropdown-item" href="/community">Community</a>
              <a class="dropdown-item" href="/stats">Statistics</a>
              <a class="dropdown-item" href="/publications">Publications</a>
              <a class="dropdown-item" href="/code_of_conduct">Code of conduct</a>
              <a class="dropdown-item" href="/join">Join nf-core</a>
            </div>
          </li>
        </ul>
        <hr class="d-md-none">
        <ul class="navbar-nav d-md-none">
          <li class="nav-item p-1">
            <a class="nav-link" target="_blank" href="https://nf-co.re/join/slack">Chat on Slack</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" target="_blank" href="https://groups.google.com/forum/#!forum/nf-core">Join the email list</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" target="_blank" href="https://twitter.com/nf_core">Follow on twitter</a>
          </li>
          <li class="nav-item p-1 mb-3">
            <a class="nav-link" target="_blank" href="https://github.com/nf-core">See nf-core on GitHub</a>
          </li>
        </ul>
        <a class="d-none d-lg-block btn btn-success" style="position:absolute; right: 1rem;" href="/join">
          Join nf-core
        </a>
      </div>
    </nav>

<?php

if(isset($title) and $title): ?>

    <div class="mainpage">
      <div class="mainpage-heading">
        <div class="container">
          <?php
          if(isset($md_github_url) and $md_github_url){
            echo '<a href="'.$md_github_url.'" class="edit-md-btn btn btn-sm btn-outline-light float-right d-none d-md-inline-block ml-2 mt-4" title="Edit this page on GitHub" data-toggle="tooltip" data-delay=\'{ "show": 500, "hide": 0 }\'><i class="fas fa-pencil-alt"></i> Edit</a>';
          }
          if(isset($header_btn_url) && isset($header_btn_text)){
            echo '<a href="'.$header_btn_url.'" class="btn btn-sm btn-outline-light float-right d-none d-md-inline-block mt-4">'.$header_btn_text.'</a>';
          }
          ?>
          <h1 class="display-3"><?php echo $title; ?></h1>
          <?php if($subtitle){ echo '<p class="lead">'.$subtitle.'</p>'; } ?>
          <?php if(isset($header_html)){ echo $header_html; } ?>
        </div>
      </div>

      <div class="triangle triangle-down"></div>

      <?php if(!isset($mainpage_container) or $mainpage_container): ?> <div class="container main-content pt-5"> <?php endif; ?>

<?php endif;
if( isset($markdown_fn) and $markdown_fn){
  // Print the parsed HTML
  if( !isset($no_print_content) or !$no_print_content ){
    echo $content;
  }
}
