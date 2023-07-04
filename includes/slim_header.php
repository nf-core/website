<?php

require_once dirname(__FILE__) . '/functions.php';

// Find the latest commit hash to prevent caching assets
$git_sha = trim(shell_exec('cd ' . dirname(__FILE__) . ' && git rev-parse --short=7 HEAD'));
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
    require_once 'parse_md.php';
    $parsed_out = parse_md($markdown_fn);
    $content = $parsed_out['content'];
    $meta = $parsed_out['meta'];
    $title = $parsed_out['title'];
    $subtitle = $parsed_out['subtitle'];
}
?>
<!doctype html>
<html lang="en">

<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="color-scheme" content="light dark">
    <!-- FontAwesome -->
    <script src="https://kit.fontawesome.com/38356a05cc.js" crossorigin="anonymous"></script>
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-68098153-2"></script>
    <!-- Other JS -->
    <script src="/assets/lib/jquery.min.js"></script>
    <script src="/assets/lib/bootstrap.bundle.min.js"></script>
    <script src="/assets/lib/highlight.min.js"></script>
    <script src="/assets/lib/groovy.min.js"></script>
    <!-- Page-specific CSS and JS -->
    <?php
    if (isset($import_moment) && $import_moment): ?>
        <script src="/assets/lib/moment.js"></script>
        <script src="/assets/lib/moment-timezone-with-data-10-year-range.min.js"></script>
    <?php endif;
    if (isset($import_chartjs) && $import_chartjs): ?>
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
    if (isset($import_schema_launcher) && $import_schema_launcher): ?>
        <script src="/assets/lib/moment.js"></script>
        <script src="/assets/lib/showdown.min.js"></script>
        <script src="/assets/js/nf-core-schema-launcher.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif;
    if (isset($aws) && $aws): ?>
        <link rel="stylesheet" href="/assets/lib/dataTables.bootstrap5.min.css">
        <script src="/assets/lib/aws-sdk.min.js"></script>
        <script src="/assets/lib/jquery.dataTables.min.js"></script>
        <script src="/assets/lib/dataTables.bootstrap5.min.js"></script>
        <script src="/assets/js/aws-s3-explorer.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif;
    if (isset($import_schema_builder) && $import_schema_builder): ?>
        <link href="/assets/lib/jquery-ui.min.css" rel="stylesheet">
        <script src="/assets/lib/jquery-ui.min.js"></script>
        <script src="/assets/lib/moment.js"></script>
        <script src="/assets/lib/showdown.min.js"></script>
        <script src="/assets/js/nf-core-schema-builder.js?c=<?php echo $git_sha; ?>"></script>
    <?php endif;
    if (isset($import_typeform) && $import_typeform): ?>
        <script src="https://embed.typeform.com/next/embed.js"></script>

    <?php endif;
    if (isset($title) && $title == 'Publications'): ?>
        <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
        <script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
    <?php endif;
    if (isset($youtube_embed)): ?>
        <script src="/assets/js/transcripts.js"></script>
    <?php endif;
    ?>

    <!-- Custom nf-core CSS and JS -->
    <link href="/assets/css/nf-core-<?php echo $theme; ?>.css?c=<?php echo $git_sha; ?>" rel="stylesheet" id="theme-stylesheet">
    <script src="/assets/js/nf-core.js?c=<?php echo $git_sha; ?>"></script>
</head>

<body tabindex="0">
    <?php
    if (isset($title) and $title): ?>
        <div class="mainpage">

            <?php if (
                !isset($mainpage_container) or $mainpage_container
            ): ?> <div class="container main-content pt-5"> <?php endif; ?>

            <?php endif;
    if (isset($markdown_fn) and $markdown_fn) {
        // Print the parsed HTML
        if (!isset($no_print_content) or !$no_print_content) {
            echo $content;
        }
    }

