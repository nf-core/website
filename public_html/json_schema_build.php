<?php

// Only keep cached schema for 24 hours
define("MAX_JSON_BUILD_CACHE_AGE", 60*60*24);

// Build URL for this page
if(isset($_SERVER['HTTPS']) && $_SERVER['HTTPS'] === 'on') $self_url = "https://";
else $self_url = "http://";
$self_url .= $_SERVER['HTTP_HOST'].$_SERVER['REQUEST_URI'];

// Holder for cache that we may load later
$schema_cache = false;
$schema = false;
$expires_timestamp = false;

function return_json($response){
    header('Content-type: application/json');
    echo json_encode($response, JSON_PRETTY_PRINT);
    exit;
}

//
// Cache house keeping
//

// Check that the cache directory exists, create it if not
$cache_dir = dirname(dirname(__FILE__)).'/api_cache/json_builder';
if (!file_exists($cache_dir)) {
    mkdir($cache_dir, 0777, true);
}

// Loop through files and delete any that are too old
$schema_cache_files = glob($cache_dir.'/*.json');
foreach($schema_cache_files as $fn) {
    $fn_parts = explode('_', basename($fn));
    if(count($fn_parts) == 2 && is_numeric($fn_parts[0])){
        $fn_expires = $fn_parts[0] + MAX_JSON_BUILD_CACHE_AGE;
        if(time() > $fn_expires){
            unlink($fn);
        }
    }
}

//
// POST request - we've been sent a JSON Schema
//
if(isset($_POST['post_content']) && $_POST['post_content'] == 'json_schema'){

    // Was a cache ID supplied?
    if(isset($_POST['cache_id'])){
        validate_cache_id($_POST['cache_id']);
        $cache_id = $_POST['cache_id'];
    } else {
        // Build a string for the filename with the timestamp and random string
        $cache_id = time().'_'.substr(md5(rand()), 0, 12);
    }
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';

    // Build a dict with the schema and 'status' => 'waiting_for_user'
    $schema_cache = array(
        'version' => $_POST['version'],
        'schema' => $_POST['schema'],
        'status' => $_POST['status']
    );
    $schema_cache_json = json_encode($schema_cache, JSON_PRETTY_PRINT)."\n";

    // Write to JSON file
    file_put_contents($cache_fn, $schema_cache_json);

    // Return with URL to check status cache if this came from the API (nf-core tools)
    if(isset($_POST['api']) && $_POST['api'] == 'true'){
        return_json(array(
            'status' => 'recieved',
            'web_url' => $self_url.'?id='.$cache_id,
            'api_url' => $self_url.'?id='.$cache_id.'&api=true'
        ));
    }

    // Not API, it came from the page form - redirect to the web url
    else {
        header('Location: '.$self_url.'?id='.$cache_id);
        exit;
    }
}

// GET request - polling for the results of a Schema builder
elseif(isset($_GET['id'])){
    $cache_id = $_GET['id'];
    validate_cache_id($cache_id);
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';

    $schema_cache = json_decode(file_get_contents($cache_fn), true);
    $schema = json_decode($schema_cache['schema'], true);

    // API check response
    if(isset($_GET['api']) && $_GET['api'] = 'true'){
        // Return just 'waiting_for_user' if flag still set
        if($schema_cache['status'] == 'waiting_for_user'){
            return_json(array('status' => 'waiting_for_user'));
        }
        // Presumably has been saved, so just return everything
        else {
            return_json($schema_cache);
        }
    }
}

function validate_cache_id($cache_id){
    // Check that ID looks valid
    $id_parts = explode('_', $cache_id);
    if(count($id_parts) != 2 || !is_numeric($id_parts[0])){
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache ID looks wrong: '.$cache_id
        ));
    }

    // Check that timestamp isn't too old
    global $expires_timestamp;
    $expires_timestamp = $id_parts[0] + MAX_JSON_BUILD_CACHE_AGE;
    if(time() > $expires_timestamp){
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache is too old (max '.(MAX_JSON_BUILD_CACHE_AGE/60/60).' hours): '.date('r', $id_parts[0])
        ));
    }

    // Check if temp file exists, return error if not
    global $cache_dir;
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';
    if (!file_exists($cache_fn)) {
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache not found: '.$cache_id
        ));
    }
}

// Got this far without printing JSON - build web GUI
$title = 'Parameter schema';
$subtitle = 'Customise a JSON Schema for your pipeline parameters';
if($schema_cache) $import_schema_builder = true;
include('../includes/header.php');

if(!$schema_cache){ ?>

<p class="mt-5">Typically this page is launched automatically and prefilled by the <code>nf-core schema build</code> command
(see the <a href="/tools">Tools</a> page for more information), however you can also paste in a JSON Schema below.</p>

<form method="post" action="">
    <input type="hidden" name="post_content" value="json_schema">
    <input type="hidden" name="version" value="web_input">
    <div class="form-group">
        <label for="schema_input">Paste your JSON Schema:</label>
        <textarea name="schema" id="schema_input" class="form-control text-monospace small" rows=10></textarea>
    </div>
    <button type="submit" class="btn btn-primary">Submit</button>
</form>

<?php } else { ?>

    <p class="lead">Schema cache ID: <code id="schema_cache_id"><?php echo $cache_id; ?></code> <small class="cache_expires_at" style="display:none;">(expires <span><?php echo $expires_timestamp; ?></span>)</small></p>
    <div class="alert alert-info small">
        <strong>Need something with more power?</strong>
        Try <a href="https://jsondraft.com/" target="_blank">jsondraft.com</a>,
        <a href="https://json-schema-editor.tangramjs.com/editor.html#/" target="_blank">json-schema-editor.tangramjs.com</a>
        or <a href="https://json-schema.org/implementations.html#editors" target="_blank">other JSON Schema editors</a>.
    </div>

    <div class="schema-builder-header sticky-top">
        <div class="row align-items-center">
            <div class="col-sm-auto">
                <button class="btn btn-outline-secondary add-param-btn"><i class="fas fa-plus-square"></i> Parameter</button>
            </div>
            <div class="col">
                <button class="btn btn-block btn-light schema-panel-btn" data-target="#schema-builder">nf-core parameter JSON Schema builder</button>
            </div>
            <div class="col-sm-auto">
                <button class="btn btn-primary schema-panel-btn" data-target="#schema-finished">Finished</button>
            </div>
        </div>
    </div>
    <div class="schema-panel" id="schema-builder"></div>
    <div class="schema-panel" id="schema-finished" style="display:none;">
        <div class="card mb-5">
            <div class="card-body">
                <h4 id="schema-send-status">Ok, that's it - done!</h4>
                <p>If you still have <code>nf-core schema build</code> running in the background,
                    it should now update with your new schema. You can close this window.</p>
                <p>If you didn't use <code>nf-core schema build</code> or it has already exited,
                    copy the schema below and paste it in to your pipeline's <code>nextflow_schema.json</code> file.</p>
                <button class="btn btn-block btn-outline-secondary copy-schema-btn">
                    <i class="far fa-copy mr-2"></i> Copy pipeline schema
                </button>
            </div>
        </div>
        <textarea id="json_schema" class="form-control text-monospace disabled" disabled rows="30"><?php echo json_encode($schema, JSON_PRETTY_PRINT); ?></textarea>
    </div>

    <!-- Params schema settings modal -->
    <div class="modal fade" id="settings_modal" tabindex="-1" role="dialog" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <span class="modal-title h4 text-monospace">params.<span></span></span>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                      <span aria-hidden="true">&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <p class="text-muted text-italic no_special_settings">No special settings available</p>
                    <div class="form-group settings_enum_group">
                        <label for="settings_enum">Enumerated values</label>
                        <input type="text" class="form-control" id="settings_enum" placeholder="value_1|value_2|value_3">
                        <small class="form-text text-muted">
                            Input values must be one of these. Separate with the pipe (<code>|</code>) character.
                            <a href="https://json-schema.org/understanding-json-schema/reference/generic.html#enumerated-values" target="_blank">
                                <i class="fas fa-question-circle"></i>
                            </a>
                        </small>
                    </div>
                    <div class="form-group settings_pattern_group">
                        <label for="settings_pattern">Pattern</label>
                        <input type="text" class="form-control" id="settings_pattern" placeholder="^[A-Za-z_][A-Za-z0-9_]*$">
                        <small class="form-text text-muted">
                            Regular expression to validate the input against
                            <a href="https://json-schema.org/understanding-json-schema/reference/string.html#regular-expressions" target="_blank">
                                <i class="fas fa-question-circle"></i>
                            </a>
                        </small>
                    </div>
                    <div class="settings_minmax_group">
                        <div class="row">
                            <div class="col-sm-6">
                                <div class="form-group mb-0">
                                    <label for="settings_minimum">Minimum</label>
                                    <input type="number" class="form-control" id="settings_minimum">
                                </div>
                            </div>
                            <div class="col-sm-6">
                                <div class="form-group mb-0">
                                    <label for="settings_maximum">Maximum</label>
                                    <input type="number" class="form-control" id="settings_maximum">
                                </div>
                            </div>
                        </div>
                        <small class="form-text text-muted mb-1">
                            Number that value must be less than / greater than or equal to
                            <a href="https://json-schema.org/understanding-json-schema/reference/numeric.html#range" target="_blank">
                                <i class="fas fa-question-circle"></i>
                            </a>
                        </small>
                    </div>
                </div>
                <div class="modal-footer row">
                    <div class="col-auto">
                        <button type="button" class="btn btn-outline-danger" data-dismiss="modal">Delete parameter</button>
                    </div>
                    <div class="col text-right">
                        <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
                        <button type="button" class="btn btn-primary" data-dismiss="modal" id="settings_save">Save changes</button>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <?php

    // TODO - write status and buttons to save schema when finished

}

include('../includes/footer.php');
