<?php
$title = 'Error 404';
$subtitle = 'Page not found';
$request_url = 'that page';
if($_SERVER['REQUEST_URI'] != '/404'){
    $request_url = '<code>https://'.$_SERVER['HTTP_HOST'].$_SERVER['REQUEST_URI'].'</code>';
}
include('../includes/header.php');
?>
<h1>Sorry, that page could not be found.</h1>
<p>It looks like <?php echo $request_url; ?> doesn't exist. </p>
<p>Please let us know by <a href="https://github.com/nf-core/nf-co.re" target="_blank">creating an issue on GitHub</a>.</p>

<?php include('../includes/footer.php'); ?>
