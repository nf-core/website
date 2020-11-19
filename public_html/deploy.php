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
  shell_exec("cd /home/nfcore/nf-co.re && git fetch && git reset --hard origin/master >> /home/nfcore/update.log 2>&1 &");

  // Pull the new version of the tools repo
  shell_exec("cd /home/nfcore/nf-co.re/includes/nf-core/tools && git fetch && git reset --hard origin/master >> /home/nfcore/update.log 2>&1 &");
  // Get the latest version number and write to a file
  shell_exec("cd /home/nfcore/nf-co.re/includes/nf-core/tools && git describe --tags --abbrev=0 > /home/nfcore/nf-co.re/includes/nf-core/tools_version.txt");
  // Print all version tags to a file
  shell_exec("cd /home/nfcore/nf-co.re/includes/nf-core/tools && git tag > /home/nfcore/nf-co.re/includes/nf-core/tools_all_versions.txt");
  // Pull the new version of the tools repo (code documentation)
  shell_exec("cd /home/nfcore/tools-docs && git fetch && git reset --hard origin/api-doc >> /home/nfcore/update.log 2>&1 &");

  die("\ndeploy done " . date("Y-m-d h:i:s") . "\n\n");

}
