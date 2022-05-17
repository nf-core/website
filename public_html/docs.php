<?php
error_reporting(E_ALL);
ini_set('display_errors', 1);

$docs_md_base = '../markdown/';

$path_parts = explode('/', $_GET['path']);
$path_parts = array_filter($path_parts); # Clear any empty array elements
$md_fn = $_GET['path'];
if (substr($md_fn, 0, 6) == 'docs/') {
    $md_fn = substr($md_fn, 6);
}
if (substr($md_fn, -3) !== '.md') {
    $md_fn .= '.md';
}
echo '<br><br><br><br><br><br>' . $md_fn;

# Render docs sub-page
if (file_exists($docs_md_base . $md_fn)) {
    $markdown_fn = $docs_md_base . $md_fn;
    $md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . $md_fn;
    $section = trim($path_parts[1]);
    include '../includes/documentation.php';
    exit();
    # Unrecognised URL - trigger 404
} elseif ($md_fn != 'index.php' && $md_fn != '') {
    header('HTTP/1.1 404 Not Found');
    include '404.php';
    die();
}

# Docs index page
$title = 'Documentation';
$subtitle = 'The nf-core documentation';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/markdown/docs';
include '../includes/header.php';
?>

hello world


<?php include '../includes/footer.php';
