<?php
// Find the latest commit hash to prevent caching assets
$git_sha = trim(shell_exec("cd ".dirname(__FILE__)." && git rev-parse --short=7 HEAD"));
if(strlen($git_sha) != 7){
  $git_sha = '';
}

?><!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>nf-core</title>
    <meta name="description" content="A collection of high quality Nextflow pipelines">
    <meta name="author" content="Phil Ewels">
    <link rel="shortcut icon" href="assets/img/logo/nf-core-logo-square.png" type="image/png" />
    <link href="assets/css/bootstrap.min.css" rel="stylesheet">
    <link href="assets/css/code_highlighting/github.css" rel="stylesheet" >
    <link href="https://use.fontawesome.com/releases/v5.0.11/css/all.css" rel="stylesheet" integrity="sha384-p2jx59pefphTFIpeqCcISO9MdVfIm4pNnsL08A6v5vaQc4owkQqxMV8kg4Yvhaw/" crossorigin="anonymous">
    <link href="assets/css/nf-core.css?c=<?php echo $git_sha; ?>" rel="stylesheet">
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-68098153-2"></script>
    <script>window.dataLayer = window.dataLayer || []; function gtag(){dataLayer.push(arguments);}  gtag('js', new Date()); gtag('config', 'UA-68098153-2'); </script>
  </head>
  <body>

    <nav class="navbar fixed-top navbar-expand-md navbar-light site-nav">
      <a class="navbar-brand d-md-none" href="/"><img height="25px" src="assets/img/logo/nf-core-logo.svg" /></a>
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
            <a class="nav-link" href="/usage_docs">Usage</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/developer_docs">Developers</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/tools">Tools</a>
          </li>
          <li class="nav-item p-1">
            <a class="nav-link" href="/about">About</a>
          </li>
        </ul>
        <hr class="d-md-none">
        <ul class="navbar-nav d-md-none">
          <li class="nav-item p-1">
            <a class="nav-link" target="_blank" href="https://gitter.im/nf-core/Lobby">Chat on Gitter</a>
          </li>
          <li class="nav-item p-1 mb-3">
            <a class="nav-link" target="_blank" href="https://github.com/nf-core">See nf-core on GitHub</a>
          </li>
        </ul>
        <div class="d-none d-md-block" style="position:absolute; right: 1rem;">
          <a class="nav-link d-inline-block px-2" target="_blank" href="https://gitter.im/nf-core/Lobby" data-toggle="tooltip" title="Chat on Gitter"><img height="25px" src="assets/img/gitter.svg" /></a>
          <a class="nav-link d-inline-block px-2" target="_blank" href="https://github.com/nf-core" data-toggle="tooltip" title="See nf-core on GitHub"><img height="25px" src="assets/img/github.svg" /></a>
        </div>
      </div>
    </nav>

<?php if(isset($title) and $title): ?>

    <div class="mainpage">
      <div class="mainpage-heading">
        <div class="container">
          <h1 class="display-3"><?php echo $title; ?></h1>
          <?php if($subtitle): ?>
          <p class="lead"><?php echo $subtitle; ?></p>
          <?php endif; ?>
        </div>
      </div>

      <div class="triangle triangle-down"></div>

      <div class="container main-content">

<?php endif;

// Convert Markdown to HTML if a filename is given
if( isset($markdown_fn) and $markdown_fn){
  // Markdown parsing libraries
  require_once('../includes/libraries/parsedown/Parsedown.php');
  require_once('../includes/libraries/parsedown-extra/ParsedownExtra.php');

  // Load the docs markdown and convert to html
  $md = file_get_contents($markdown_fn);
  $pd = new ParsedownExtra();
  $content = $pd->text($md);

  // Print the parsed HTML
  if( !isset($no_print_content) or !$no_print_content ){
    echo $content;
  }
}
?>
