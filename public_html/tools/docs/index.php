<?php
$title = 'nf-core/tools - Documentation';
$subtitle = 'The nf-core companion tool, to help with common tasks.';
include '../../../includes/header.php';
?>
<p class="lead">This page shows the available documentation for the "back-end" of the nf-core/tools package.</p>
<p>It's primarily for two audiences: pipeline developers wanting to know more about how the pipeline lint tests function,
    and Python developers using the nf-core/tools Python library in their own packages.</p>

<p><blockquote>To read documentation aimed at command-line users of nf-core/tools, please visit <a href="../">the main tools page</a>.</blockquote></p>

<p>Documentation is built using Sphinx and is available below for each released version of tools:</p>

<?php
$dirs = array_filter(glob('[0-9]*'), 'is_dir');
usort($dirs, 'version_compare');
echo '<ul>';
foreach (array_merge(['latest', 'dev'], array_reverse($dirs)) as $dir) {
    echo '<li><a href="' . $dir . '">' . $dir . "</a></li>\n";
}
echo '</ul>';

include '../../../includes/footer.php';


?>
