<?php // Common functionality for JSON Schema - build and launch

require_once 'functions.php';

// Only keep cached schema for 2 weeks
define('MAX_JSON_BUILD_CACHE_AGE', 60 * 60 * 24 * 14);

// Get URL for this page
$self_url = get_self_url();

// Holder for cache that we may load later
$cache = false;
$schema = false;
$expires_timestamp = false;

//
// Cache house keeping
//

// Check that the cache directory exists, create it if not
if (!file_exists($cache_dir)) {
    mkdir($cache_dir, 0777, true);
}

// Loop through files and delete any that are too old
$cache_files = glob($cache_dir . '/*.json');
foreach ($cache_files as $fn) {
    $fn_parts = explode('_', basename($fn));
    if (count($fn_parts) == 2 && is_numeric($fn_parts[0])) {
        $fn_expires = $fn_parts[0] + MAX_JSON_BUILD_CACHE_AGE;
        if (time() > $fn_expires) {
            unlink($fn);
        }
    }
}

//
// POST request - we've been sent a JSON Schema
//
if (isset($_POST['post_content']) && $_POST['post_content'] == $post_content_type) {
    // Was a cache ID supplied?
    if (isset($_POST['cache_id'])) {
        validate_cache_id_api($_POST['cache_id']);
        $cache_id = $_POST['cache_id'];
        $cache_id = basename($cache_id);
    } else {
        // Build a string for the filename with the timestamp and random string
        $cache_id = time() . '_' . substr(md5(rand()), 0, 12);
    }
    $cache_fn = $cache_dir . '/' . $cache_id . '.json';

    // Build a dict with the schema and 'status' => 'waiting_for_user'
    $cache = [];
    foreach ($post_keys as $k) {
        $cache[$k] = $_POST[$k];
    }
    $cache_json = json_encode($cache, JSON_PRETTY_PRINT) . "\n";

    // Write to JSON file
    file_put_contents($cache_fn, $cache_json);

    // Return with URL to check status cache if this came from the API (nf-core tools)
    if (isset($_POST['api']) && $_POST['api'] == 'true') {
        return_json([
            'status' => 'recieved', // DO NOT FIX THIS TYPO. nf-core/tools will break.
            'web_url' => $self_url . '?id=' . $cache_id,
            'api_url' => $self_url . '?id=' . $cache_id . '&api=true',
        ]);
    }

    // Not API, it came from the page form - redirect to the web url
    else {
        header('Location: '.strtok($_SERVER['REQUEST_URI'], '?') .'?id=' . $cache_id);
        exit();
    }
}

// GET request - polling for the results of a Schema builder
elseif (isset($_GET['id']) && !isset($_POST['post_content'])) {
    $cache_id = $_GET['id'];
    $cache_id = basename($cache_id);
    validate_cache_id_api($cache_id);
    $cache_fn = $cache_dir . '/' . $cache_id . '.json';

    $cache = json_decode(file_get_contents($cache_fn), true);

    // Decode JSON objects
    foreach (['schema', 'nxf_flags', 'input_params'] as $k) {
        if (isset($cache[$k])) {
            $parsed_json = json_decode($cache[$k], true);
            if ($parsed_json !== null) {
                $cache[$k] = $parsed_json;
            }
        }
    }

    // API check response
    if (isset($_GET['api']) && ($_GET['api'] = 'true')) {
        // Return just 'waiting_for_user' if flag still set
        if ($cache['status'] == 'waiting_for_user') {
            return_json(['status' => 'waiting_for_user']);
        }
        // Presumably has been saved, so just return everything
        else {
            return_json($cache);
        }
    }
}

function validate_cache_id_api($cache_id) {
    $check = validate_cache_id($cache_id);
    if ($check['status'] == 'error') {
        return_json($check);
    }
}

function validate_cache_id($cache_id) {
    // Check that ID looks valid
    $id_parts = explode('_', $cache_id);
    if (count($id_parts) != 2 || !is_numeric($id_parts[0])) {
        return [
            'status' => 'error',
            'message' => 'JSON Build cache ID looks wrong: ' . $cache_id,
        ];
    }

    // Check that timestamp isn't too old
    global $expires_timestamp;
    $expires_timestamp = $id_parts[0] + MAX_JSON_BUILD_CACHE_AGE;
    if (time() > $expires_timestamp) {
        return [
            'status' => 'error',
            'message' =>
                'JSON Build cache is too old (max ' .
                MAX_JSON_BUILD_CACHE_AGE / 60 / 60 / 24 .
                ' days): ' .
                date('r', $id_parts[0]),
        ];
    }

    // Check if temp file exists, return error if not
    global $cache_dir;
    $cache_fn = $cache_dir . '/' . $cache_id . '.json';
    if (!file_exists($cache_fn)) {
        return [
            'status' => 'error',
            'message' => 'JSON Build cache not found: ' . $cache_id,
        ];
    }
    return [
        'status' => 'success',
        'message' => 'JSON Build cache looks valid: ' . $cache_id,
    ];
}
