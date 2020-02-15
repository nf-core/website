<?php
$title = 'Publications';
$subtitle = 'Publications about the nf-core project and its pipelines';
$markdown_fn = '../markdown/publications.md';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/markdown/publications.md';
$no_print_content = true;
include('../includes/header.php');

$altmetric_pattern = '/<!-- pub-stats (\S+) -->/';
$altmetric_html = '
<div class="pub-stats-wrapper">
    <div data-doi="${1}" data-badge-popover="bottom" data-badge-type="donut" data-hide-no-mentions="true" class="altmetric-embed"></div>
</div>
<div class="pub-stats-wrapper">
    <span data-doi="${1}" class="__dimensions_badge_embed__" data-hide-zero-citations="true" data-style="small_circle" data-legend="hover-bottom"></span>
</div>';
$content = preg_replace($altmetric_pattern, $altmetric_html, $content);
echo $content;

echo '
<!-- Dimensions.ai -->
<script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
<!-- Altmetric -->
<script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
';

include('../includes/footer.php');
