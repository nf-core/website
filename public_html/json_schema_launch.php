<?php

$cache_dir = dirname(dirname(__FILE__)).'/api_cache/json_launch';
$post_content_type = 'json_schema_launcher';
$post_keys = ['version', 'schema', 'nxf_flags', 'input_params', 'status'];
require_once('../includes/json_schema.php');

// Got this far without printing JSON - build web GUI
$title = 'Launch pipeline';
$subtitle = 'Configure workflow parameters for a pipeline run';
if($cache) $import_schema_launcher = true;
$mainpage_container = false;
include('../includes/header.php');
?>
<div class="container">

<?php if(!$cache){ ?>

<h3>Launch a pipeline</h3>

<p>You can run <code>nf-core launch</code> to submit any pipeline schema to this page and set the parameters required for launch.</p>

<?php } else { ?>

    <p class="lead">Params cache ID: <code id="params_cache_id"><?php echo $cache_id; ?></code> <small class="cache_expires_at" style="display:none;">(expires <span><?php echo $expires_timestamp; ?></span>)</small></p>
</div>
<div class="container-fluid main-content">
    <h1>Not yet built</h1>
    <h3>Set status to saved</h3>
    <p id="schema-send-status"></p>
    <button class="btn btn-primary launcher-panel-btn" data-target="#params-finished">Set as saved</button>
    <h3>Cache results</h3>
    <pre><?php print_r($cache); ?></pre>
    <h3>Schema</h3>
    <textarea id="json_schema" class="form-control text-monospace disabled" disabled rows="30"><?php echo json_encode($schema, JSON_PRETTY_PRINT); ?></textarea>
</div> <!-- .container-fluid -->

<?php } // if $cache

include('../includes/footer.php');
