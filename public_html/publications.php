<?php
$title = 'Publications';
$subtitle = 'Publications about the nf-core project and its pipelines';
$markdown_fn = '../markdown/publications.md';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/markdown/publications.md';
$no_print_content = true;
include('../includes/header.php');

$altmetric_insert = '<!-- nf-core-altmetric -->';
$altmetric_donut = '<div style="float:right; background-color:#ffffff; border-radius:5px; padding:15px 15px 0 15px; margin:-15px -15px 0 0;" data-badge-details="right" data-badge-type="medium-donut" data-doi="10.1038/s41587-020-0439-x" data-hide-no-mentions="true" class="altmetric-embed"></div>';
$content = str_replace($altmetric_insert, $altmetric_donut, $content);
echo $content;

echo '
<!-- Altmetric -->
<script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
';

include('../includes/footer.php');
