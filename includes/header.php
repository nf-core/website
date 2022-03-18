<?php

require_once(dirname(__FILE__) . '/functions.php');

// Find the latest commit hash to prevent caching assets
$git_sha = trim(shell_exec("cd " . dirname(__FILE__) . " && git rev-parse --short=7 HEAD"));
if (strlen($git_sha) != 7) {
  $git_sha = '';
}

// Theme switcher cookie
$theme = 'auto';
if (isset($_COOKIE['nfcoretheme']) && in_array($_COOKIE['nfcoretheme'], ['auto', 'light', 'dark'])) {
  $theme = $_COOKIE['nfcoretheme'];
}

// Convert Markdown to HTML if a filename is given
if (isset($markdown_fn) and $markdown_fn) {
  require_once('parse_md.php');
  $parsed_out = parse_md($markdown_fn);
  $content = $parsed_out["content"];
  $meta = $parsed_out["meta"];
  $title = $parsed_out["title"];
  $subtitle = $parsed_out["subtitle"];
}

// Page title
$page_title = 'nf-core';
if (isset($title) && strlen($title) > 0) {
  $page_title = preg_replace('/^nf-core\//', '', strip_tags($title)) . ' &raquo; nf-core';
}
// Page meta description
$page_meta = 'A collection of high quality Nextflow pipelines';
if (isset($subtitle) && strlen($subtitle) > 0) {
  $page_meta = strip_tags($subtitle);
}

?>
<!doctype html>
<html lang="en">

<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
  <title><?php echo $page_title; ?></title>
  <meta name="description" content="<?php echo $page_meta; ?>">
  <meta name="author" content="Phil Ewels">
  <meta name="color-scheme" content="light dark">
  <link rel="shortcut icon" href="/assets/img/logo/nf-core-logo-square.png" type="image/png" />
  <link rel="alternate" type="application/rss+xml" title="nf-core: Events" href="/events/rss" />
  <!-- FontAwesome -->
  <script src="https://kit.fontawesome.com/471b59d3f8.js"></script>
  <!-- Global site tag (gtag.js) - Google Analytics -->
  <script async src="https://www.googletagmanager.com/gtag/js?id=UA-68098153-2"></script>
  <!-- Other JS -->
  <script src="/assets/lib/jquery.min.js"></script>
  <script src="/assets/lib/bootstrap.bundle.min.js"></script>
  <script src="/assets/lib/highlight.min.js"></script>
  <script src="/assets/lib/groovy.min.js"></script>
  <!-- Page-specific CSS and JS -->
  <?php if (isset($import_moment) && $import_moment) : ?>
    <script src="/assets/lib/moment.js"></script>
    <script src="/assets/lib/moment-timezone-with-data-10-year-range.min.js"></script>
  <?php endif;
  if (isset($import_leaflet) && $import_leaflet) : ?>
    <link href="/assets/lib/leaflet.css" rel="stylesheet">
    <script src="/assets/lib/leaflet.js"></script>
  <?php endif;
  if (isset($import_chartjs) && $import_chartjs) : ?>
    <script src="/assets/lib/moment.js"></script>
    <script src="/assets/lib/Chart.min.js"></script>
    <link href="/assets/lib/Chart.min.css" rel="stylesheet">
    <script src="/assets/lib/hammer.min.js"></script>
    <script src="/assets/lib/chartjs-plugin-zoom.min.js"></script>
    <!-- NO CURRENT VERSION ON NPM -->
    <script src="/assets/js/canvas2svg.js"></script><!-- NO CURRENT VERSION ON NPM -->
    <script src="/assets/lib/FileSaver.min.js"></script>
    <script src="/assets/lib/jquery.tablesorter.min.js"></script>
  <?php endif;
  if (isset($import_schema_launcher) && $import_schema_launcher) : ?>
    <script src="/assets/lib/moment.js"></script>
    <script src="/assets/lib/showdown.min.js"></script>
    <script src="/assets/js/nf-core-schema-launcher.js?c=<?php echo $git_sha; ?>"></script>
  <?php endif;
  if (isset($aws) && $aws) : ?>
    <link rel="stylesheet" href="/assets/lib/dataTables.bootstrap5.min.css">
    <script src="/assets/lib/aws-sdk.min.js"></script>
    <script src="/assets/lib/jquery.dataTables.min.js"></script>
    <script src="/assets/lib/dataTables.bootstrap5.min.js"></script>
    <script src="/assets/js/aws-s3-explorer.js?c=<?php echo $git_sha; ?>"></script>
  <?php endif;
  if (isset($import_schema_builder) && $import_schema_builder) : ?>
    <link href="/assets/lib/jquery-ui.min.css" rel="stylesheet">
    <script src="/assets/lib/jquery-ui.min.js"></script>
    <script src="/assets/lib/moment.js"></script>
    <script src="/assets/lib/showdown.min.js"></script>
    <script src="/assets/js/nf-core-schema-builder.js?c=<?php echo $git_sha; ?>"></script>
  <?php endif;
  if (isset($import_typeform) && $import_typeform) : ?>
    <script src="https://embed.typeform.com/next/embed.js"></script>

  <?php endif;
  if (isset($title) && $title == 'Publications') : ?>
    <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
    <script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
  <?php endif;
  if (isset($youtube_embed)) : ?>
    <script src="https://www.youtube.com/iframe_api"></script>
    <script src="/assets/js/transcripts.js"></script>
  <?php endif; ?>

  <!-- Custom nf-core CSS and JS -->
  <link href="/assets/css/nf-core-<?php echo $theme; ?>.css?c=<?php echo $git_sha; ?>" rel="stylesheet" id="theme-stylesheet">
  <script src="/assets/js/nf-core.js?c=<?php echo $git_sha; ?>"></script>
  <script>
    window.dataLayer = window.dataLayer || [];

    function gtag() {
      dataLayer.push(arguments);
    }
    gtag('js', new Date());
    gtag('config', 'UA-68098153-2');

    <?php if (isset($youtube_embed)) : ?>
      youtube_embed = true;
    <?php endif; ?>
  </script>
</head>

<body tabindex="0">
  <nav class="navbar fixed-top navbar-expand-md navbar-light bg-light shadow-sm site-nav d-print-none">
    <div class="container-fluid">
      <a class="navbar-brand d-md-none" href="/">
        <img height="25px" src="/assets/img/logo/nf-core-logo.svg" class="hide-dark">
        <img height="25px" src="/assets/img/logo/nf-core-logo-darkbg.svg" class="hide-light">
      </a>
      <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarCollapse" aria-controls="navbarCollapse" aria-expanded="false" aria-label="Toggle navigation">
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
            <a class="nav-link" href="/modules">Modules</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/tools">Tools</a>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link dropdown-toggle" href="/usage/introduction" role="button" data-bs-toggle="dropdown">Usage</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/usage/introduction">Getting started</a>
              <a class="dropdown-item" href="/usage/installation">Installation</a>
              <a class="dropdown-item" href="/usage/configuration">Pipeline configuration</a>
              <a class="dropdown-item" href="/usage/offline">Running offline</a>
              <a class="dropdown-item" href="/usage/usage_tutorials">Usage tutorials</a>
              <a class="dropdown-item" href="/usage/reference_genomes">Reference genomes</a>
              <a class="dropdown-item" href="/usage/data_management">Data Management</a>
              <a class="dropdown-item" href="/usage/troubleshooting">Troubleshooting</a>
              <a class="dropdown-item" href="/usage/nextflow">Nextflow resources</a>
            </div>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link dropdown-toggle" href="/developers/guidelines" role="button" data-bs-toggle="dropdown">Developers</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/developers/guidelines">Guidelines</a>
              <a class="dropdown-item" href="/developers/adding_pipelines">Adding a new pipeline</a>
              <a class="dropdown-item" href="/developers/modules">Modules</a>
              <a class="dropdown-item" href="/developers/release_checklist">Release checklist</a>
              <a class="dropdown-item" href="/tools/docs/latest/pipeline_lint_tests/">Lint error codes</a>
              <a class="dropdown-item" href="/developers/sync">Template synchronisation</a>
              <a class="dropdown-item" href="/developers/developer_tutorials">Developer tutorials</a>
              <a class="dropdown-item" href="/developers/editor_plugins">Code editor plugins</a>
              <a class="dropdown-item" href="/developers/design_guidelines">Graphic design guidelines</a>
            </div>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/events"><?php if ($curr_event and $curr_event['ongoing']) {
                                                  echo '<i class="fad fa-circle text-danger me-1"></i>';
                                                } ?>Events</a>
          </li>
          <li class="nav-item p-1 dropdown">
            <a class="nav-link dropdown-toggle" href="/about" role="button" data-bs-toggle="dropdown">About</a>
            <div class="dropdown-menu">
              <a class="dropdown-item" href="/about">About nf-core</a>
              <a class="dropdown-item" href="/community">Community</a>
              <a class="dropdown-item" href="/stats">Statistics</a>
              <a class="dropdown-item" href="/publications">Publications</a>
              <a class="dropdown-item" href="/mentorships">Mentorships</a>
              <a class="dropdown-item" href="/code_of_conduct">Code of conduct</a>
              <a class="dropdown-item" href="/join">Join nf-core</a>
            </div>
          </li>
        </ul>
        <hr class="d-md-none">
        <a class="d-md-none btn d-block btn-success mb-3" href="/join">
          Join nf-core
        </a>
        <a class="d-none d-lg-block btn btn-success" style="position:absolute; right: 1rem;" href="/join">
          Join nf-core
        </a>
      </div>
    </div>
  </nav>

  <?php

  if (isset($title) and $title) : ?>

    <div class="mainpage">
      <div class="mainpage-heading triangle-down">
        <div class="container">
          <?php
          if (isset($md_github_url) and $md_github_url) {
            echo '<a href="' . $md_github_url . '" class="edit-md-btn btn btn-sm btn-outline-light float-end d-none d-md-inline-block ms-2 mt-4 d-print-none" title="Edit this page on GitHub" data-bs-toggle="tooltip" data-bs-delay=\'{ "show": 500, "hide": 0 }\'><i class="fas fa-pencil-alt"></i> Edit</a>';
          }
          if (isset($header_btn_url) && isset($header_btn_text)) {
            echo '<a href="' . $header_btn_url . '" class="btn btn-sm btn-outline-light float-end d-none d-md-inline-block mt-4">' . $header_btn_text . '</a>';
          }
          ?>
          <h1 class="display-2"><?php echo $title; ?></h1>
          <?php if ($subtitle) {
            echo '<p class="lead">' . $subtitle . '</p>';
          } ?>
          <?php if (isset($header_html)) {
            echo $header_html;
          } ?>
        </div>
      </div>

      <?php if (!isset($mainpage_container) or $mainpage_container) : ?> <div class="container main-content pt-5"> <?php endif; ?>

      <?php endif;
    if (isset($markdown_fn) and $markdown_fn) {
      // Print the parsed HTML
      if (!isset($no_print_content) or !$no_print_content) {
        echo $content;
      }
    }
