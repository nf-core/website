<?php
$title = 'Publications';
$subtitle = 'Publications about the sanger-tol project and its pipelines';
$markdown_fn = '../markdown/publications.md';
$md_github_url = 'https://github.com/sanger-tol/pipelines-website/blob/main/markdown/publications.md';

$altmetric_pattern = '/<!-- pub-stats (\S+) -->/';
$altmetric_html = '
<div class="pub-stats-wrapper ms-2 mt-3 me-1">
    <div data-doi="${1}" data-badge-popover="bottom" data-badge-type="donut" data-hide-no-mentions="true" class="altmetric-embed"></div>
</div>
<div class="pub-stats-wrapper ms-2 mt-3 me-1">
    <span data-doi="${1}" class="__dimensions_badge_embed__" data-hide-zero-citations="true" data-style="small_circle" data-legend="hover-bottom"></span>
</div>';
$html_content_replace = [$altmetric_pattern, $altmetric_html];
include '../includes/header.php';
include '../includes/footer.php';
