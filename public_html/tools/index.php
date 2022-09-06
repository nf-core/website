<?php
$title = 'Tools';
$subtitle = 'The nf-core companion tool, to help with common tasks.';

$markdown_fn = dirname(__FILE__) . '/../../markdown/tools/README.md'; // Local cache
$md_github_url = ' https://github.com/nf-core/tools/blob/master/README.md'; // For the footer
$md_github_raw_url = 'https://raw.githubusercontent.com/nf-core/tools/dev/README.md'; // For rendering the page
$src_url_prepend = 'https://raw.githubusercontent.com/nf-core/tools/master/'; // For rich-codex generated images

// Fetch readme if cache is not found or more than 24 hours old
if (!file_exists($markdown_fn) || filemtime($markdown_fn) < time() - 60 * 60 * 24) {
    // Download the readme and cache
    // Build directories if needed
    if (!is_dir(dirname($markdown_fn))) {
        mkdir(dirname($markdown_fn), 0777, true);
    }
    $md_contents = file_get_contents($md_github_raw_url);
    if ($md_contents) {
        file_put_contents($markdown_fn, $md_contents);
    }
}

$md_trim_before = '## Table of contents';
include '../../includes/header.php';
include '../../includes/footer.php';
?>
