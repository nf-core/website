<?php
$title = 'nf-core/tools - Documentation';
$subtitle = 'A companion tool is available for nf-core to help with common tasks.';
include('../../../includes/header.php');

echo '<p>Documentation for each version of nf-core/tools:</p><ul>';
$dirs = array_filter(glob('*'), 'is_dir');
foreach($dirs as $dir){
    echo '<li><a href="'.$dir.'">'.$dir."</a></li>\n";
}
echo '</ul>';

include('../../../includes/footer.php');
?>
