<?php
if ( $_POST['payload'] ) {

  // Get the hook secret
  $config = parse_ini_file("../config.ini");

  // Check that the hashed signature looks right
  list($algo, $hash) = explode('=', $_SERVER['HTTP_X_HUB_SIGNATURE'], 2) + array('', '');
  $rawPost = file_get_contents('php://input');
  if ($hash !== hash_hmac($algo, $rawPost, $config['secret'])){
    die('GitHub signature looks wrong.');
  }

  // Pull the new version of the website repo
  shell_exec("cd /home/nfcore/nf-co.re && git fetch && git reset --hard origin/master");
  // Pull the new version of the tools repo
  shell_exec("cd /home/nfcore/nf-co.re/includes/nf-core/tools && git fetch && git reset --hard origin/master");
  // Update the JSON file describing all the pipeline versions
  shell_exec("php /home/nfcore/nf-co.re/update_pipeline_details.php");
  die("done " . mktime());

}
?>
