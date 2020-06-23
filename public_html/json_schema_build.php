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
$mainpage_container = false;
include('../includes/header.php');
?>
<div class="container">

<p class="mt-5">nf-core pipelines have a file in their root directories called <code>nextflow_schema.json</code> which
describes the input parameters that the pipeline accepts.
This page helps pipeline authors to build their pipeline schema file by using a graphical interface.</p>

<div class="row">
    <div class="col"><hr></div>
    <div class="col-auto"><a class="text-muted" data-toggle="collapse" href="#page_help" role="button"><small>click here to read more</small></a></div>
    <div class="col"><hr></div>
</div>
<div class="collapse" id="page_help">
    <h5>nf-core schema</h5>

    <p>nf-core schema files use the <a href="https://json-schema.org/" target="_blank">JSON Schema</a> standard,
    with a couple of extra assumptions:</p>
    <ul>
        <li>The top-level schema is an <code>object</code>, where each of the <code>properties</code> are a pipeline parameter or group of parameters</li>
        <li>Only groups one-<code>object</code> deep are allowed (groups or no groups are fine, but no nested groups)</li>
        <li>Before use in validation, all groups are flattened. As such, no duplicate keys are allowed across groups.</li>
    </ul>
    <p>We also use a couple of extra JSON keys:</p>
    <ul>
        <li><code>help_text</code>, a longer description providing more in-depth help. Typically <code>description</code> is just a few words long and the longer help text is shown when a user requests it.</li>
        <li><code>hidden: True</code>, which tells tools to ignore this <code>param</code> in interfaces and help text etc.</li>
        <li><code>fa_icon</code>, a <a href="https://fontawesome.com/" target="_blank">fontawesome.com</a> icon for use in web interfaces (eg: <code>&lt;i class="fas fa-flask"&gt;&lt;/i&gt;</code> - <i class="fas fa-flask"></i> )</li>
    </ul>

    <h5>Tips:</h5>
    <ul>
        <li>Click the <i class="fas fa-cog"></i> icon on the right to access more settings.</li>
        <li>Click and drag the <i class="fas fa-grip-vertical"></i> icon on the left to re-order parameters and groups.</li>
        <li>The <i class="fas fa-icons"></i> allows you to set a custom icon for the parameter or group.</li>
        <li>The <i class="fas fa-book help_text_icon"></i> icon shows whether help text has been written. To add, click on it.</li>
        <li>Be a power user with keyboard shortcuts!
            <ul class="small">
                <li>Use <code class="border shadow-sm">Enter</code> and <code class="border shadow-sm">Shift</code>+<code class="border shadow-sm">Enter</code> to go up and down.</li>
                <li>Use <code class="border shadow-sm">Tab</code> and <code class="border shadow-sm">Shift</code>+<code class="border shadow-sm">Tab</code> to go right and left.</li>
                <li><code class="border shadow-sm">Space</code> toggles checkboxes and opens dropdown boxes.</li>
                <li><code class="border shadow-sm">ctrl</code>+<code class="border shadow-sm">shift</code>+<code class="border shadow-sm">,</code> opens the settings panel.</li>
                <li><code class="border shadow-sm">ctrl</code>+<code class="border shadow-sm">shift</code>+<code class="border shadow-sm">&uarr;</code> moves the row up.</li>
                <li><code class="border shadow-sm">ctrl</code>+<code class="border shadow-sm">shift</code>+<code class="border shadow-sm">&darr;</code> moves the row down.</li>
            </ul>
        </li>
    </ul>

    <div class="row">
        <div class="col"><hr></div>
        <div class="col-auto"><a class="text-muted" data-toggle="collapse" href="#page_help" role="button"><small>hide read-more text</small></a></div>
        <div class="col"><hr></div>
    </div>
</div>

<?php if(!$schema_cache){ ?>

<h3>Load Schema</h3>

<p>If you previously ran <code>nf-core schema build</code> and forgot to save, you can resume editing by entering the build ID below:</p>
<form method="get" action="" class="form-inline mb-2">
    <input type="text" class="form-control mr-2" name="id" placeholder="Build ID">
    <button type="submit" class="btn btn-primary">Load</button>
</form>

<p>Note that if you did save and want to continue editing, just run <code>nf-core schema build</code> in your workflow again.</p>

<h3>New Schema</h3>

<p>Typically this page is launched automatically and prefilled by the <code>nf-core schema build</code> command
(see the <a href="/tools">Tools</a> page for more information), however you can also paste in a JSON Schema below.</p>

<form method="post" action="">
    <input type="hidden" name="post_content" value="json_schema">
    <input type="hidden" name="version" value="web_input">
    <div class="form-group">
        <label for="schema_input">Paste your JSON Schema:</label>
        <textarea name="schema" id="schema_input" class="form-control text-monospace small" rows=10>
{
    "$schema": "http://json-schema.org/draft-07/schema",
    "$id": "https://raw.githubusercontent.com/YOUR_PIPELINE/master/nextflow_schema.json",
    "title": "Nextflow pipeline parameters",
    "description": "This pipeline uses Nextflow and processes some kind of data. The JSON Schema was built using the nf-core pipeline schema builder.",
    "type": "object",
    "properties": {
        "some_parameter": {
            "type": "string"
        }
    }
}
        </textarea>
    </div>
    <button type="submit" class="btn btn-primary mb-3">Submit</button>
</form>

<?php } else { ?>

    <p class="lead">Schema cache ID: <code id="schema_cache_id"><?php echo $cache_id; ?></code> <small class="cache_expires_at" style="display:none;">(expires <span><?php echo $expires_timestamp; ?></span>)</small></p>
</div>
<div class="container-fluid main-content">

    <div class="schema-builder-header sticky-top">
        <div class="row align-items-center">
            <div class="col-8 col-lg-4">
                <button class="btn btn-outline-secondary add-param-btn"><i class="fas fa-plus-square mr-1"></i> Parameter</button>
                <button class="btn btn-outline-secondary add-group-btn"><i class="fas fa-folder-plus mr-1"></i> Group</button>
            </div>
            <div class="col-lg-4 d-none d-lg-block">
                <button class="btn btn-block btn-light schema-panel-btn" data-target="#schema-builder">
                    nf-core schema builder
                </button>
            </div>
            <div class="col-4 text-right">
                <button class="btn btn-primary schema-panel-btn" data-target="#schema-finished"><i class="fas fa-check-square mr-1"></i> Finished</button>
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
                <button class="btn btn-block btn-outline-secondary copy-schema-btn m-auto">
                    <i class="far fa-copy mr-2"></i> Copy pipeline schema
                </button>
            </div>
        </div>
    </div>

    <h3>Pipeline JSON Schema</h3>
    <p>This is the schema for your pipeline. As you change values in the form above, it will update. When you are finished, click <em>Finished</em> in the top toolbar.</p>
    <textarea id="json_schema" class="form-control text-monospace disabled" disabled rows="30"><?php echo json_encode($schema, JSON_PRETTY_PRINT); ?></textarea>

    <!-- Params schema settings modal -->
    <div class="modal fade" id="settings_modal" tabindex="-1" role="dialog" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <span class="modal-title h4 text-monospace"></span>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                      <span aria-hidden="true">&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <p class="text-muted settings_nothing_special">No special settings available.</p>
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
                <div class="modal-footer">
                    <div class="col-auto pl-0">
                        <button type="button" class="btn btn-outline-danger" data-dismiss="modal" id="settings_delete"><i class="fas fa-trash-alt mr-1"></i> <span>Delete parameter</span></button>
                    </div>
                    <div class="col text-right pr-0">
                        <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
                        <button type="button" class="btn btn-primary" data-dismiss="modal" id="settings_save">Save changes</button>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- Help text modal -->
    <div class="modal fade" id="help_text_modal" tabindex="-1" role="dialog" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <span class="modal-title h4 text-monospace"></span>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                      <span aria-hidden="true">&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <label for="help_text_input">
                        Help text is used for generating documentation and is shown on demand for the command-line and web launch tools.
                    </label>
                    <div class="card">
                        <div class="card-header">
                            <ul class="nav nav-tabs card-header-tabs">
                                <li class="nav-item"><a class="nav-link active" data-toggle="tab" href="#tab-helptext">Write</a></li>
                                <li class="nav-item"><a class="nav-link" data-toggle="tab" href="#tab-helptext-preview">Preview</a></li>
                            </ul>
                        </div>
                        <div class="card-body tab-content">
                            <div class="tab-pane fade show active" id="tab-helptext" role="tabpanel">
                                <div class="form-group help_text_modal_group mb-0">
                                    <textarea class="form-control" id="help_text_input" rows="5"></textarea>
                                    <small class="form-text text-muted mt-2">
                                        <i class="fab fa-markdown"></i> <a href="https://www.markdownguide.org/cheat-sheet/" target="_blank">Markdown</a> supported,
                                        but this will be shown raw on the command line so please keep it simple.
                                    </small>
                                </div>
                            </div>
                            <div class="tab-pane fade" id="tab-helptext-preview" role="tabpanel">
                                <p>Command-line:</p>
                                <pre><span class="helptext-cli-preview-title"></span>
<span class="helptext-preview-description"></span>

<span class="helptext-preview-helptext text-muted"></span></pre>
                                <p>Website:</p>
                                <div class="card helptext-html-preview">
                                    <div class="card-body">
                                        <h4 class="mt-0 pt-0 helptext-web-preview-title"></h4>
                                        <p class="lead helptext-preview-description"></p>
                                        <div class="helptext-preview-helptext"></div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="modal-footer">
                    <div class="col text-right">
                        <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
                        <button type="button" class="btn btn-primary" data-dismiss="modal" id="help_text_save">Save</button>
                    </div>
                </div>
            </div>
        </div>
    </div>

</div> <!-- .container-fluid -->

<?php } // if $schema_cache

include('../includes/footer.php');
