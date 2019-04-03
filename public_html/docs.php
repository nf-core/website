<?php
$md_base = dirname(dirname(__file__))."/markdown/";
$md_fn = $_GET['md'];
# Add on the .md file extension
if(substr($md_fn, -3) !== '.md'){
    $md_fn .= '.md';
}
if(file_exists($md_base.$md_fn)){
    $markdown_fn = $md_base.$md_fn;
}
# .htaccess rewriting truncates these directories
else if(file_exists($md_base.'usage/'.$md_fn)){
    $markdown_fn = $md_base.'usage/'.$md_fn;
} else if(file_exists($md_base.'developers/'.$md_fn)){
    $markdown_fn = $md_base.'developers/'.$md_fn;
} else {
    header('Location: /404');
}
include('../includes/header.php');
include('../includes/footer.php');
?>
