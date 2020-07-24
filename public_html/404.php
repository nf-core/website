<?php
header('HTTP/1.1 404 Not Found');

$title = 'Error 404';
$subtitle = 'Page not found';
$request_url = 'that page';
if($_SERVER['REQUEST_URI'] != '/404'){
  $protocol = isset($_SERVER['HTTPS']) && !empty($_SERVER['HTTPS']) ? 'https://' : 'http://';
  $request_url = '<code>'.$protocol.$_SERVER['HTTP_HOST'].$_SERVER['REQUEST_URI'].'</code>';
}
include('../includes/header.php');
?>
<h1>Sorry, that page could not be found.</h1>
<p>It looks like <?php echo $request_url; ?> doesn't exist. </p>
<p>Please let us know by <a href="https://github.com/nf-core/nf-co.re" target="_blank">creating an issue on GitHub</a>.</p>

<?php
if(isset($suggestion_404_url)){
  echo '<p>You might be trying to get to: <a href="'.$protocol.$_SERVER['HTTP_HOST'].$suggestion_404_url.'">'.$protocol.$_SERVER['HTTP_HOST'].$suggestion_404_url.'</a>';
}

include('../includes/footer.php'); ?>
