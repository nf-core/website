<?php
$title = 'Publications';
$subtitle = 'Publications about the nf-core project and its pipelines';
$markdown_fn = '../markdown/publications.md';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/markdown/publications.md';
$no_print_content = true;
$no_auto_toc = true;
include('../includes/header.php');

$altmetric_pattern = '/<!-- pub-stats (\S+) -->/';
$altmetric_html = '
<div class="pub-stats-wrapper ml-2 mt-3 mr-1">
    <div data-doi="${1}" data-badge-popover="bottom" data-badge-type="donut" data-hide-no-mentions="true" class="altmetric-embed"></div>
</div>
<div class="pub-stats-wrapper ml-2 mt-3 mr-1">
    <span data-doi="${1}" class="__dimensions_badge_embed__" data-hide-zero-citations="true" data-style="small_circle" data-legend="hover-bottom"></span>
</div>';
$content = preg_replace($altmetric_pattern, $altmetric_html, $content);

echo '<div class="row flex-wrap-reverse flex-lg-wrap"><div class="col-12 col-lg-9">';

echo '<div class="rendered-markdown publication-page-content">' . $content . '</div>';

echo '</div>'; # end of the content col
echo '<div class="col-12 col-lg-3 pl-2"><div class="side-sub-subnav sticky-top">';
# ToC
$toc = '<nav class="toc">';
$toc .= generate_toc($content);
# Back to top link
$toc .= '<p class="small text-right"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
$toc .= '</nav>';
echo $toc;
echo '</div></div>'; # end of the sidebar col
echo '</div>'; # end of the row

include('../includes/footer.php');
