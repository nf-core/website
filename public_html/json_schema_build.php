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
    $fn_parts = explode('_', $fn);
    if(count($fn_parts) == 2 && is_numeric($fn_parts[0])){
        if($fn_parts[0] < (time() - MAX_JSON_BUILD_CACHE_AGE)){
            unlink($fn);
        }
    }
}

//
// POST request - we've been sent a JSON Schema
//
if(isset($_POST['post_content']) && $_POST['post_content'] == 'json_schema'){

    // Build a string for the filename with the timestamp and random string
    $cache_id = time().'_'.substr(md5(rand()), 0, 12);
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';

    // Build a dict with the schema and 'status' => 'waiting_for_user'
    $schema_cache = array(
        'version' => $_POST['version'],
        'schema' => $_POST['schema'],
        'status' => 'waiting_for_user'
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

// TODO - build code to handle schema submitted from buider tool

// GET request - polling for the results of a Schema builder
elseif(isset($_GET['id'])){
    $cache_id = $_GET['id'];

    // Check that ID looks valid
    $id_parts = explode('_', $cache_id);
    if(count($id_parts) != 2 || !is_numeric($id_parts[0])){
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache ID looks wrong: '.$cache_id
        ));
    }

    // Check that timestamp isn't too old
    $expires_timestamp = $id_parts[0] + MAX_JSON_BUILD_CACHE_AGE;
    if(time() > $expires_timestamp){
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache is too old (max '.(MAX_JSON_BUILD_CACHE_AGE/60/60).' hours): '.date('r', $id_parts[0])
        ));
    }

    // Check if temp file exists, return error if not
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';
    if (!file_exists($cache_fn)) {
        return_json(array(
            'status' => 'error',
            'message' => 'JSON Build cache not found: '.$cache_id
        ));
    }

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

// Got this far without printing JSON - build web GUI
$title = 'JSON Schema Builder';
$subtitle = 'Customise a JSON Schema for your pipeline';
$import_schema_builder = true;
include('../includes/header.php');

if(!$schema_cache){ ?>

<p class="mt-5">Typically this page is launched automatically and prefilled by the <code>nf-core schema build</code> command
(see the <a href="/tools">Tools</a> page for more information), however you can also paste in a JSON Schema below.</p>

<form method="post" action="">
    <input type="hidden" name="post_content" value="json_schema">
    <input type="hidden" name="version" value="web">
    <div class="form-group">
        <label for="schema_input">Paste your JSON Schema:</label>
        <textarea name="schema" id="schema_input" class="form-control text-monospace small" rows=10></textarea>
    </div>
    <button type="submit" class="btn btn-primary">Submit</button>
</form>

<?php } else { ?>

    <p class="lead">Schema cache ID: <code><?php echo $cache_id; ?></code> <small class="cache_expires_at" style="display:none;">(expires <span><?php echo $expires_timestamp; ?></span>)</small></p>
    <div class="alert alert-info small">
        <strong>Need something with more power?</strong>
        Try <a href="https://jsondraft.com/" target="_blank">jsondraft.com</a>,
        <a href="https://json-schema-editor.tangramjs.com/editor.html#/" target="_blank">json-schema-editor.tangramjs.com</a>
        or <a href="https://json-schema.org/implementations.html#editors" target="_blank">other JSON Schema editors</a>.
    </div>

    <div class="schema-builder"></div>

    <pre id="json_schema"><?php echo json_encode($schema, JSON_PRETTY_PRINT); ?></pre>

    <?php

    // TODO - write status and buttons to save schema when finished

}

include('../includes/footer.php');
