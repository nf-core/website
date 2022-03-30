<?php
if ($_POST['payload']) {
    // Get the hook secret
    $config = parse_ini_file('../config.ini');

    // Check that the hashed signature looks right
    [$algo, $hash] = explode('=', $_SERVER['HTTP_X_HUB_SIGNATURE'], 2) + ['', ''];
    $rawPost = file_get_contents('php://input');
    if ($hash !== hash_hmac($algo, $rawPost, $config['secret'])) {
        die('GitHub signature looks wrong.');
    }

    // Update the JSON file describing all the pipeline versions
    shell_exec('php ' . dirname(dirname(__FILE__)) . '/update_pipeline_details.php >> /home/nfcore/update.log 2>&1 &');

    die("\ndeploy_pipelines done " . date('Y-m-d h:i:s') . "\n\n");
}
