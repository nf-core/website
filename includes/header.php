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
if( isset($markdown_fn) and $markdown_fn){
  // Markdown parsing libraries
  require_once(dirname(__FILE__).'/libraries/parsedown/Parsedown.php');
  require_once(dirname(__FILE__).'/libraries/parsedown-extra/ParsedownExtra.php');
  require_once(dirname(__FILE__).'/libraries/Spyc.php');

  // Load the docs markdown
  $md_full = file_get_contents($markdown_fn);
  if ($md_full === false) {
    header('HTTP/1.1 404 Not Found');
    include('404.php');
    die();
  }

  // Get the meta
  $meta = [];
  $md = $md_full;
  if(substr($md_full,0,3) == '---'){
    $md_parts = explode('---', $md_full, 3);
    if(count($md_parts) == 3){
      $meta = spyc_load($md_parts[1]);
      $md = $md_parts[2];
      if(isset($meta['title'])){
        $title = $meta['title'];
      }
      if(isset($meta['subtitle'])){
        $subtitle = $meta['subtitle'];
      }
    }
  }

  // Trim off any content if requested
  if(isset($md_trim_before) && $md_trim_before){
    // Only trim if the string exists
    if(stripos($md, $md_trim_before)){
      $md = stristr($md, $md_trim_before);
    }
  }
  if(isset($md_trim_after) && $md_trim_after){
    if(stripos($md, $md_trim_after)){
      $md = stristr($md, $md_trim_after);
    }
  }

  // Find and replace markdown content if requested
  if(isset($md_content_replace)){
    $md = str_replace($md_content_replace[0], $md_content_replace[1], $md);
  }

  // Format Nextflow code blocks as Groovy
  $md = preg_replace('/```nextflow/i', '```groovy', $md);

  // Convert to HTML
  $pd = new ParsedownExtra();
  $content = $pd->text($md);

  // Highlight any search terms if we have them
  if(isset($_GET['q']) && strlen($_GET['q'])){
    $content = preg_replace("/(".$_GET['q'].")/i", "<mark>$1</mark>", $content);
  }

  // Automatically add HTML IDs to headers
  // Add ID attributes to headers
  $hids = Array();
  $content = preg_replace_callback(
    '~<h([1234])>(.*?)</h([1234])>~Ui', // Ungreedy by default, case insensitive
    function ($matches) {
      global $hids;
      $id_match = strip_tags($matches[2]);
      $id_match = strtolower( preg_replace('/[^\w\-\.]/', '', str_replace(' ', '-', $id_match)));
      $id_match = str_replace('---', '-', $id_match);
      $hid = $id_match;
      $i = 1;
      while(in_array($hid, $hids)){
        $hid = $id_match.'-'.$i;
        $i += 1;
      }
      $hids[] = $hid;
      return '<h'.$matches[1].' id="'.$hid.'"><a href="#'.$hid.'" class="header-link"><span class="fas fa-link"></span></a>'.$matches[2].'</h'.$matches[3].'>';
    },
    $content);

    // Prepend to src URLs if configureds and relative
    if(isset($src_url_prepend)){
      $content = preg_replace('/src="(?!https?:\/\/)([^"]+)"/i', 'src="'.$src_url_prepend.'$1"', $content);
    }
    // Prepend to href URLs if configureds and relative
    if(isset($href_url_prepend)){
      $content = preg_replace('/href="(?!https?:\/\/)(?!#)([^"]+)"/i', 'href="'.$href_url_prepend.'$1"', $content);
    }
    // Clean up href URLs if configured
    if(isset($href_url_suffix_cleanup)){
      $content = preg_replace('/href="(?!https?:\/\/)(?!#)([^"]+)'.$href_url_suffix_cleanup.'"/i', 'href="$1"', $content);
    }
    // Find and replace HTML content if requested
    if(isset($html_content_replace)){
      $content = str_replace($html_content_replace[0], $html_content_replace[1], $content);
    }

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
    <link href="/assets/css/leaflet.css" rel="stylesheet">
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
    <script src="/assets/js/leaflet.js"></script>
    <!-- Page-specific CSS and JS -->
    <?php if(isset($import_chartjs) && $import_chartjs): ?>
    <script src="/assets/js/moment.js"></script>
    <script src="/assets/js/Chart.min.js"></script>
    <script src="/assets/js/hammer.min.js"></script>
    <script src="/assets/js/chartjs-plugin-zoom.min.js"></script>
    <script src="/assets/js/canvas2svg.js"></script>
    <script src="/assets/js/FileSaver.js"></script>
    <?php endif; ?>
    <?php if(isset($import_schema_builder) && $import_schema_builder): ?>
    <link href="/assets/css/jquery-ui.min.css" rel="stylesheet">
    <link href="/assets/css/iconpicker.css" rel="stylesheet">
    <script src="/assets/js/jquery-ui.min.js"></script>
    <script src="/assets/js/moment.js"></script>
    <script src="/assets/js/iconpicker.js"></script>
    <script src="/assets/js/nf-core-schema-builder.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif; ?>
    <script src="/assets/js/jquery.tablesorter.min.js"></script>

    <!-- Custom nf-core CSS and JS -->
    <link href="/assets/css/nf-core.css?c=<?php echo $git_sha; ?>" rel="stylesheet">
    <link href="/assets/css/nf-core-<?php echo $theme; ?>.css?c=<?php echo $git_sha; ?>" rel="stylesheet" id="theme-stylesheet">
    <script src="/assets/js/nf-core.js?c=<?php echo $git_sha; ?>"></script>
    <script>window.dataLayer = window.dataLayer || []; function gtag(){dataLayer.push(arguments);}  gtag('js', new Date()); gtag('config', 'UA-68098153-2'); </script>
  </head>
  <body>

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
          <?php if(isset($header_btn_url) && isset($header_btn_text)){
            echo '<a href="'.$header_btn_url.'" class="btn btn-outline-light float-right d-none d-md-inline-block mt-4">'.$header_btn_text.'</a>';
          } ?>
          <h1 class="display-3"><?php echo $title; ?></h1>
          <?php if($subtitle){ echo '<p class="lead">'.$subtitle.'</p>'; } ?>
          <?php if(isset($header_html)){ echo $header_html; } ?>
        </div>
      </div>

      <div class="triangle triangle-down"></div>

      <?php if(!isset($mainpage_container) or $mainpage_container): ?> <div class="container main-content"> <?php endif; ?>

<?php endif;
if( isset($markdown_fn) and $markdown_fn){
  // Print the parsed HTML
  if( !isset($no_print_content) or !$no_print_content ){
    echo $content;
  }
}
