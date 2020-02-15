<?php
$title = 'Publications';
$subtitle = 'Publications about the nf-core project and its pipelines';
$markdown_fn = '../markdown/publications.md';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/markdown/publications.md';
$no_print_content = true;
include('../includes/header.php');

$altmetric_pattern = '/<!-- altmetric (\S+) -->/';
$altmetric_html = '<div class="altmetric-wrapper"><div data-badge-popover="left" data-badge-type="donut" data-doi="${1}" data-hide-no-mentions="true" class="altmetric-embed"></div></div>';
$content = preg_replace($altmetric_pattern, $altmetric_html, $content);
echo $content;

echo '
<!-- Altmetric -->
<script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
';

include('../includes/footer.php');
