<?php

require_once('../includes/functions.php');
$error_msgs = array();

// Get available pipelines / releases
$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = array();
foreach($pipelines_json->remote_workflows as $wf){
    if($wf->archived) continue;
    $releases = [];
    if(count($wf->releases) > 0){
        usort($wf->releases, 'rsort_releases');
        foreach($wf->releases as $release){
            $releases[] = $release->tag_name;
        }
    }
    $releases[] = 'dev';
    $pipelines[$wf->name] = $releases;
}

// Loading launch page for a pipeline from the website
$nxf_flag_schema = array(
    'Nextflow command-line flags' => [
        'type' => 'object',
        'description' => 'General Nextflow flags to control how the pipeline runs.',
        'help_text' => "These are not specific to the pipeline and will not be saved in any parameter file. They are just used when building the `nextflow run` launch command.",
        'properties' => [
            '-name' => [
                'type' => 'string',
                'description' => 'Unique name for this nextflow run',
                'pattern' => '^[a-zA-Z0-9-_]+$'
            ],
            '-profile' => [
                'type' => 'string',
                'description' => 'Configuration profile'
            ],
            '-work-dir' => [
                'type' => 'string',
                'description' => 'Work directory for intermediate files',
                'default' => './work',
            ],
            '-resume' => [
                'type' => 'boolean',
                'description' => 'Resume previous run, if found',
                'help_text' => "Execute the script using the cached results, useful to continue executions that was stopped by an error",
                'default' => False
            ]
        ]
    ]
);
if(isset($_GET['pipeline']) && isset($_GET['release'])){
    $error_msgs = launch_pipeline_web($_GET['pipeline'], $_GET['release']);
}
function launch_pipeline_web($pipeline, $release){
    // Check that we recognise the pipeline name
    global $pipelines;
    global $nxf_flag_schema;
    if(!array_key_exists($_GET['pipeline'], $pipelines)){
        return ["Error - Pipeline name <code>$pipeline</code> not recognised"];
    }
    // Try to fetch the nextflow_schema.json file
    $gh_pipeline_schema_fn = dirname(dirname(__FILE__))."/api_cache/json_schema/{$pipeline}/{$release}.json";
    $gh_pipeline_no_schema_fn = dirname(dirname(__FILE__))."/api_cache/json_schema/{$pipeline}/{$release}.NO_SCHEMA";
    # Build directories if needed
    if (!is_dir(dirname($gh_pipeline_schema_fn))) {
      mkdir(dirname($gh_pipeline_schema_fn), 0777, true);
    }
    // Load cache if not 'dev'
    if(file_exists($gh_pipeline_no_schema_fn) && $release != 'dev'){
        return [
            "Error - Could not find a pipeline schema for <code>$pipeline</code> - <code>$release</code>",
            "Please launch using the command line tool instead: <code>nf-core launch $pipeline -r $release</code>",
            "<!-- URL attempted: $gh_launch_schema_url -->"
        ];
    } else if(file_exists($gh_pipeline_schema_fn) && $release != 'dev'){
        $gh_launch_schema_json = file_get_contents($gh_pipeline_schema_fn);
    } else {
        $api_opts = stream_context_create([ 'http' => [ 'method' => 'GET', 'header' => [ 'User-Agent: PHP' ] ] ]);
        $gh_launch_schema_url = "https://api.github.com/repos/nf-core/{$pipeline}/contents/nextflow_schema.json?ref={$release}";
        $gh_launch_schema_json = file_get_contents($gh_launch_schema_url, false, $api_opts);
        if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
            # Remember for next time
            file_put_contents($gh_pipeline_no_schema_fn, '');
            return [
                "Error - Could not find a pipeline schema for <code>$pipeline</code> - <code>$release</code>",
                "Please launch using the command line tool instead: <code>nf-core launch $pipeline -r $release</code>",
                "<!-- URL attempted: $gh_launch_schema_url -->"
            ];
        } else {
            # Save cache
            file_put_contents($gh_pipeline_schema_fn, $gh_launch_schema_json);
        }
    }
    // Build the POST data
    $gh_launch_schema_response = json_decode($gh_launch_schema_json, true);
    $gh_launch_schema = json_decode(base64_decode($gh_launch_schema_response['content']), true);
    $gh_launch_schema['properties'] = $nxf_flag_schema + $gh_launch_schema['properties'];
    $_POST['post_content'] = 'json_schema_launcher';
    $_POST['api'] = 'false';
    $_POST['version'] = 'web_launcher';
    $_POST['status'] = 'waiting_for_user';
    $_POST['cli_launch'] = false;
    $_POST['nxf_flags'] = "{}";
    $_POST['input_params'] = "{}";
    $_POST['pipeline'] = 'nf-core/'.$pipeline;
    $_POST['revision'] = $release;
    $_POST['nextflow_cmd'] = "nextflow run $pipeline -r $release";
    $_POST['schema'] = json_encode($gh_launch_schema);
    return [];
}

// Share code to go through POST data and handle cache
$cache_dir = dirname(dirname(__FILE__)).'/api_cache/json_launch';
$post_content_type = 'json_schema_launcher';
$post_keys = ['version', 'schema', 'nxf_flags', 'input_params', 'status', 'cli_launch', 'nextflow_cmd', 'pipeline', 'revision'];
require_once('../includes/json_schema.php');

// Return to editor
if(isset($_GET['return_to_editor']) && $_GET['return_to_editor'] == 'true'){
    $cache['cli_launch'] = false;
    $cache['version'] = 'web_launcher';
    $cache['status'] = 'waiting_for_user';

    // Write to JSON file
    $cache['schema'] = json_encode($cache['schema']);
    $cache['nxf_flags'] = json_encode($cache['nxf_flags']);
    $cache['input_params'] = json_encode($cache['input_params']);
    $cache_json = json_encode($cache, JSON_PRETTY_PRINT)."\n";
    file_put_contents($cache_fn, $cache_json);
    // Redirect to web URL
    header('Location: '.$self_url.'?id='.$cache_id);
    exit;
}

// Save form output
if(isset($_POST['post_content']) && $_POST['post_content'] == "json_schema_launcher_webform"){
    $error_msgs = save_launcher_form();
}
function save_launcher_form(){
    // global vars
    global $cache_dir;
    global $cache;
    global $self_url;
    // Check cache ID
    if(!isset($_POST['cache_id'])){
        return ["No cache ID supplied"];
    }
    $id_check = validate_cache_id($_POST['cache_id']);
    if(!isset($id_check['status'])){
        return ["Problem loading cache: <pre>".$id_check.'</pre>'];
    }
    if($id_check['status'] == 'error'){
        return ["Problem loading cache: ".$id_check['message']];
    }
    // Load cache
    $cache_id = $_POST['cache_id'];
    $cache_fn = $cache_dir.'/'.$cache_id.'.json';
    if(!file_exists($cache_fn)) {
        return ["Cache file not found: <code>".$cache_fn.'</code>'];
    }
    $cache = json_decode(file_get_contents($cache_fn), true);
    if(!isset($cache['schema'])){
        return ["Cache had no schema: <code>".$cache_fn.'</code>'];
    }
    $cache['schema'] = json_decode($cache['schema'], true);
    if(!isset($cache['schema']['properties']) || count($cache['schema']['properties']) == 0){
        return ["Cache schema was empty: <code>".$cache_fn.'</code><pre>'.print_r($cache['schema'], true).'</pre>'];
    }

    // Overwrite some keys (not schema)
    $cache['version'] = 'web_launcher';
    $cache['status'] = 'launch_params_complete';
    $cache['nxf_flags'] = array();
    $cache['input_params'] = array();

    // Loop through POST vars and set params
    foreach ($_POST as $k => $v){
        if(substr($k,0,7) == 'params_'){
            $cache['input_params'][substr($k,7)] = $v;
        }
        if(substr($k,0,9) == 'nxf_flag_'){
            $cache['nxf_flags'][substr($k,9)] = $v;
        }
    }
    // Write to JSON file
    $cache['schema'] = json_encode($cache['schema']);
    $cache['nxf_flags'] = json_encode($cache['nxf_flags']);
    $cache['input_params'] = json_encode($cache['input_params']);
    $cache_json = json_encode($cache, JSON_PRETTY_PRINT)."\n";
    file_put_contents($cache_fn, $cache_json);
    // Redirect to web URL
    header('Location: '.$self_url.'?id='.$cache_id);
    exit;
}

// Markdown parsing libraries
require_once('../includes/libraries/parsedown/Parsedown.php');
require_once('../includes/libraries/parsedown-extra/ParsedownExtra.php');
$pd = new ParsedownExtra();

function parse_md($text){
    global $pd;
    // Remove whitespace on lines that are only whitespace
    $text = preg_replace('/^\s*$/m', '', $text);
    // Remove global text indentation
    $indents = array();
    foreach(explode("\n",$text) as $l){
        if(strlen($l) > 0){
            $indents[] = strlen($l) - strlen(ltrim($l));
        }
    }
    if(min($indents) > 0){
        $text = preg_replace('/^\s{'.min($indents).'}/m', '', $text);
    }
    return $pd->text($text);
}

function build_form_param($param_id, $param, $is_required){

    global $cache;

    $dash_param_id = substr($param_id, 0, 1) == '-' ? $param_id : '--'.$param_id;
    $form_param_name = 'params_'.$param_id;
    if(substr($param_id,0,1) == '-'){
        $form_param_name = 'nxf_flag_'.$param_id;
    }

    // Hidden
    $hide_class = '';
    if(isset($param['hidden']) && (strtolower($param['hidden']) == 'true' || $param['hidden'] === true)){
        $hide_class = 'is_hidden';
    }

    // Icon
    $fa_icon = '';
    if(isset($param['fa_icon'])){
        $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw mr-3"></i>';
    }

    // Description
    $description = '';
    if(isset($param['description'])){
        $description = '<small class="form-text">'.parse_md($param['description']).'</small>';
    }

    // Help text
    $help_text_btn = '';
    $help_text = '';
    if(isset($param['help_text']) && strlen(trim($param['help_text'])) > 0){
        $help_text_btn = '<div class="input-group-append" title="Show help text" data-toggle="tooltip">
            <button class="btn input-group-btn" type="button" data-toggle="collapse" href="#help-text-'.$param_id.'" aria-expanded="false">
                <i class="fas fa-question-circle"></i>
            </button>
        </div>';
        $help_text = '<div class="collapse" id="help-text-'.$param_id.'">
            <div class="card card-body small text-muted launcher-help-text">
                '.parse_md($param['help_text']).'
            </div>
        </div>';
    }

    // Schema default value
    $placeholder = '';
    $value = '';
    if(isset($param['default'])){
        $placeholder = 'placeholder="'.$param['default'].'"';
        $value = $param['default'];
    }
    // Supplied value
    if(isset($cache['input_params'][$param_id])){
        $value = $cache['input_params'][$param_id];
    }

    // Required
    $required = '';
    $required_asterisk = '';
    $validation_text = '';
    if($is_required){
        $required = 'required';
        $required_asterisk = '<sup class="text-danger ml-2" title="Required" data-toggle="tooltip">*</sup>';
        $validation_text = '<div class="invalid-feedback">This parameter is required</div>';
    }

    // Text, number, integer, range input
    $input_type = 'text';
    $step = '';
    $minimum = '';
    $maximum = '';
    $pattern = '';
    if($param['type'] == 'number' || $param['type'] == 'integer'){
        $input_type = 'number';
    }
    if($param['type'] == 'range'){
        $input_type = 'range';
    }
    if($param['type'] == 'integer'){
        $step = 'step="1"';
        $pattern = 'pattern="\d+"';
        $validation_text = '<div class="invalid-feedback">Must be an integer</div>';
    }
    if(array_key_exists('minimum', $param) && strlen($param['minimum']) > 0){
        $minimum = 'min="'.$param['minimum'].'"';
    }
    if(array_key_exists('maximum', $param) && strlen($param['maximum']) > 0){
        $maximum = 'max="'.$param['maximum'].'"';
    }
    if(array_key_exists('pattern', $param) && strlen($param['pattern']) > 0){
        $pattern = 'pattern="'.$param['pattern'].'"';
        $validation_text = '<div class="invalid-feedback">Must match pattern <code>'.$param['pattern'].'</code></div>';
    }
    $input_el = '<input type="'.$input_type.'" '.$step.' '.$minimum.' '.$maximum.' '.$pattern.' class="form-control text-monospace" id="'.$form_param_name.'" name="'.$form_param_name.'" '.$placeholder.' value="'.$value.'" '.$required.'>';

    // Boolean input
    if($param['type'] == 'boolean'){
        $input_el = '
        <div class="form-control pl-4">
            <div class="form-check form-check-inline mr-4">
                <input '.($value === true || strtolower($value) == 'true' ? 'checked' : '').' class="form-check-input" type="radio" id="'.$form_param_name.'_true" name="'.$form_param_name.'" '.$required.' value="true">
                <label class="form-check-label" for="'.$form_param_name.'_true">True</label>
            </div>
            <div class="form-check form-check-inline">
                <input '.($value === false || strtolower($value) == 'false' ? 'checked' : '').' class="form-check-input" type="radio" id="'.$form_param_name.'_false" name="'.$form_param_name.'" '.$required.' value="false">
                <label class="form-check-label" for="'.$form_param_name.'_false">False</label>
            </div>
        </div>';
    }

    // enum input
    if(array_key_exists('enum', $param) && count($param['enum']) > 0){
        $input_el = '<select class="custom-select" id="'.$form_param_name.'" name="'.$form_param_name.'" '.$required.'>';
        foreach($param['enum'] as $option){
            $input_el .= '<option '.($value == $option ? 'selected' : '').' value="'.$option.'">'.$option.'</option>';
        }
        $input_el .= '</select>';
    }

    // Build HTML
    return '
    <div class="form-group param-form-group '.$hide_class.'" id="'.$param_id.'_group">
        <div class="input-group">
            <div class="input-group-prepend">
                <label class="input-group-text text-monospace" for="'.$form_param_name.'">'.$fa_icon.$dash_param_id.$required_asterisk.'</label>
            </div>
            '.$input_el.$help_text_btn.$validation_text.'
        </div>
        '.$description.$help_text.'
    </div>';
}

// Got this far without printing JSON - build web GUI
$title = 'Launch pipeline';
$subtitle = 'Configure workflow parameters for a pipeline run.';
$import_schema_launcher = true;
include('../includes/header.php');

if(count($error_msgs) > 0){
    echo '<div class="alert alert-danger">'.implode('<br>', $error_msgs).'</div>';
}

if(!$cache){ ?>

<p class="lead mt-5">This tool shows the available parameters for a pipeline in form for you to fill in.
    It typically works in combination with the <a href="/tools"><code>nf-core</code> helper package</a>.
</p>

<h2>Launch a pipeline from the web</h2>

<p>Pick a pipeline and release below to show the launch form for that pipeline.
When you click <em>Finished</em>, your inputs will be saved and you'll be shown the commands
to use to launch the pipeline with your choices.</p>

<div class="card card-body mb-3">
    <form action="" method="get" class="row" id="launch_select_pipeline">
        <div class="col">
            <div class="form-group mb-0 mr-3">
                <div class="input-group">
                    <div class="input-group-prepend">
                        <span class="input-group-text text-monospace">nf-core launch</span>
                    </div>
                    <select class="custom-select" name="pipeline" id="launch-pipeline-name">
                        <option value="">Select a pipeline</option>
                        <option value="">--</option>
                    <?php
                    foreach($pipelines as $wf_name => $releases_json){
                        echo '<option data-releases=\''.json_encode($releases_json).'\'>'.$wf_name.'</option>';
                    }
                    ?>
                    </select>
                </div>
            </div>
        </div>
        <div class="col">
            <div class="form-group mb-0 mr-3">
                <div class="input-group">
                    <div class="input-group-prepend">
                        <span class="input-group-text text-monospace">-revision</span>
                    </div>
                    <select class="custom-select" name="release" id="launch-pipeline-release" disabled>
                        <option>Please select a pipeline</option>
                    </select>
                </div>
            </div>
        </div>
        <div class="col-auto">
            <button type="submit" class="btn btn-primary btn-launch" id="launch-pipeline-submit" disabled>
                <i class="fad fa-rocket-launch"></i> Launch
            </button>
        </div>
    </form>
</div>

<p>For more options, such as launching a custom pipeline or using a GitHub branch, please use this tool locally - see below.</p>
<p>Read more about the different nf-core pipelines on the <a href="/pipelines">Pipelines</a> page.</p>

<h2>Launch a pipeline locally</h2>

<p>You can run <code>nf-core launch</code> to submit any pipeline schema to this page and set the parameters required for launch.
This should work with any Nextflow pipeline (though the experience is best for pipelines that have a <code>nextflow_schema.json</code> file).</p>

<p>For example, to launch the <a href="/atacseq">nf-core/atacseq</a> pipeline in your current directory:</p>
<pre>nf-core launch atacseq</pre>

<p>To launch your own custom pipeline that you have locally:</p>
<pre>nf-core launch ./my_pipeline/</pre>

<p>The tool will check the pipeline's schema and create one if none exists, and then ask if you want to use this web tool or the command-line wizard:<p>
<pre>$ nf-core launch atacseq

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 1.10.dev0


INFO: This tool ignores any pipeline parameter defaults overwritten by Nextflow config files or profiles


INFO: Using local workflow: nf-core/atacseq (dev - 00d035c)

INFO: <span style="color:green;">[✓] Pipeline schema looks valid</span>

<span style="color:purple;">Would you like to enter pipeline parameters using a web-based interface or a command-line wizard?</span>

? Choose launch method (Use arrow keys)
 <span style="color:green;">&#x276f;</span> Web based
   Command line</pre>

<p>If you select <code>Web based</code>, then this web page will load with the pipeline parameters for you to fill in.
    The command-line tool will wait for you to click <em>Finished</em> and then offer to run Nextflow with the supplied parameters.</p>

<p>The command-line wizard uses the exact same procedure, but runs entirely locally with a prompt system.
    That is best for those running on offline systems, or if you are concerned about sending sensitive information over the web.</p>

<p>The nf-core website stores a cached copy of your answers for 2 weeks under a random ID.</p>

<?php } else if($cache['status'] == 'launch_params_complete') {

    $nxf_flags = ' ';
    foreach($cache['nxf_flags'] as $key => $val){
        if(!$nxf_flag_schema['Nextflow command-line flags']['properties'][$key]['default'] == $val){
            $nxf_flags .= "$key $val ";
        }
    }
    ?>

<h1>Launch parameters saved</h1>

<?php if(isset($cache['cli_launch']) && $cache['cli_launch']): ?>
<p>The <code>nf-core launch</code> command running in your terminal should have automatically detected your settings.
    Follow the prompts in your command line to launch the pipeline.</p>
<p>If the launch command has stopped running for any reason, you can still launch by following the instructions below:</p>
<?php else: ?>
<p>Your workflow parameters are ready to go! Follow the instructions below for instructions on how to launch your pipeline:</p>
<?php endif; ?>

<h3>If your system has an internet connection</h3>
<p>The easiest way to launch this workflow is by using the <code>nf-core/tools</code> helper package.</p>
<p>Once installed (<a href="https://nf-co.re/tools#installation" target="_blank">see documentation</a>),
    simply run the following command and follow the prompts:</p>
<pre>nf-core launch --id <?php echo $cache_id; ?></pre>

<h3>Launching with no internet and without nf-core/tools</h3>
<p>You can run this pipeline with just Nextflow installed by copying the JSON below to a file called <code>nf-params.json</code>:</p>
<pre><?php echo json_encode($cache['input_params']); ?></pre>

<p>Then, launch Nextflow with the following command:</p>
<pre><?php echo $cache['nextflow_cmd']; echo $nxf_flags; ?>-params-file nf-params.json</pre>

<h3>Continue editing</h3>
<p>If you would like to continue editing your workflow parameters, click the button below:</p>
<form action="" method="get">
    <input type="hidden" name="id" value="<?php echo $cache_id; ?>">
    <input type="hidden" name="return_to_editor" value="true">
    <button type="submit" class="btn btn-outline-primary">
        <i class="fas fa-pencil-alt mr-1"></i>
        Return to editor
    </button>
</form>


<?php } else {
    $pipeline_name_header = '';
    if(isset($cache['pipeline']) && strlen($cache['pipeline']) > 0 && $cache['pipeline'] != '.'){
        $pipeline_name_header = '<p class="lead">Pipeline: <code>'.$cache['pipeline'].'</code>';
        if(isset($cache['revision']) && strlen($cache['revision']) > 0){
            $pipeline_name_header .= ' (<code>'.$cache['revision'].'</code>)';
        }
        $pipeline_name_header .= '</p>';
    }
    ?>

    <div class="alert alert-info">
        <?php echo $pipeline_name_header; ?>
        <p class="lead mb-0">Launch ID: <code><?php echo $cache_id; ?></code> <small class="cache_expires_at" style="display:none;">(expires <span><?php echo $expires_timestamp; ?></span>)</small></p>
    </div>

    <p>Go through the pipeline inputs below, setting them to the values that you would like.
        When you're done, click <em>Launch</em> and your parameters will be saved.
        <?php if($cache['cli_launch']): ?>
        The <code>nf-core launch</code> command running in the background should detect your changes and give you further instructions.
        <?php endif; ?>
    </p>
    <?php if(!$cache['cli_launch']): ?>
    <p>The page shown will show a command that you can use to directly launch the workflow.
        For those running on a system with no internet connection, you can copy the parameters JSON to a file
        and use the supplied command to launch.</p>
    <?php endif; ?>

    <form id="schema_launcher_form" action="" method="post" class="needs-validation" novalidate>

        <input type="hidden" name="cache_id" value="<?php echo $cache_id; ?>">
        <input type="hidden" name="post_content" value="json_schema_launcher_webform">

        <div class="schema-gui-header sticky-top">
            <div class="row align-items-center">
                <div class="col-md-auto">
                    <div class="btn-group">
                        <button class="btn btn-outline-secondary dropdown-toggle" type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                            <i class="far fa-stream mr-1"></i> <span>Jump to section</span>
                        </button>
                        <div class="dropdown-menu" aria-labelledby="dropdownMenuButton">
                            <?php
                            foreach($cache['schema']['properties'] as $param_id => $param){
                                if($param['type'] == 'object'){
                                    $html_id = preg_replace('/[^a-z0-9-_]/', '_', preg_replace('/\s+/', '_', strtolower($param_id)));
                                    $hidden_class = 'is_hidden';
                                    foreach($cache['schema']['properties'][$param_id]['properties'] as $child_param_id => $child_param){
                                        if(!isset($child_param['hidden']) || (strtolower($child_param['hidden']) == 'false' || $child_param['hidden'] === false)){
                                            $hidden_class = '';
                                        }
                                    }
                                    $fa_icon = '';
                                    if(isset($param['fa_icon'])){
                                        $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw mr-3 text-secondary"></i>';
                                    }
                                    echo '<a class="dropdown-item '.$hidden_class.' scroll_to_link" href="#'.$html_id.'">'.$fa_icon.$param_id.'</a>';
                                }
                            }
                            ?>
                        </div>
                    </div>
                    <button class="btn btn-outline-secondary btn-show-hidden-fields" title="Parameters that do not typically need to be altered for a normal run are hidden by default" data-toggle="tooltip" data-delay='{ "show": 500, "hide": 0 }'>
                        <span class="is_not_hidden"><i class="fas fa-eye-slash mr-1"></i> Show hidden params</span>
                        <span class="is_hidden"><i class="fas fa-eye mr-1"></i> Hide hidden params</span>
                    </button>
                </div>
                <div class="col d-none d-lg-block">
                    <span id="progress_section" class="text-muted">Nextflow command-line flags</span>
                    <div class="progress" style="height: 2px;">
                        <div class="progress-bar" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"></div>
                    </div>
                </div>
                <div class="col-md-auto">
                    <button type="submit" class="btn btn-primary btn-launch" title="Save parameters and <?php if($cache['cli_launch']) echo 'return to the command line"'; else echo 'copy command to launch' ?>" data-toggle="tooltip" data-delay='{ "show": 800, "hide": 0 }'>
                        <i class="fad fa-rocket-launch mr-1"></i> Launch
                    </button>
                </div>
            </div>
        </div>

        <?php
        foreach($cache['schema']['properties'] as $param_id => $param){
            if($param['type'] == 'object'){
                $html_id = preg_replace('/[^a-z0-9-_]/', '_', preg_replace('/\s+/', '_', strtolower($param_id)));
                $hidden_class = 'is_hidden';
                $child_parameters = '';
                foreach($cache['schema']['properties'][$param_id]['properties'] as $child_param_id => $child_param){
                    $child_parameters .= build_form_param($child_param_id, $child_param, @in_array($child_param_id, $cache['schema']['properties'][$param_id]['required']));
                    if(!isset($child_param['hidden']) || (strtolower($child_param['hidden']) == 'false' || $child_param['hidden'] === false)){
                        $hidden_class = '';
                    }
                }
                $fa_icon = '';
                if(isset($param['fa_icon'])){
                    $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw mr-3"></i>';
                }
                $description = '';
                if(isset($param['description'])){
                    $description = '<p class="form-text">'.$param['description'].'</p>';
                }
                $helptext = '';
                if(isset($param['help_text'])){
                    $helptext = '<small class="form-text text-muted">'.$param['help_text'].'</small>';
                }
                if(strlen($child_parameters) > 0){
                    echo '
                    <fieldset class="'.$hidden_class.'" id="'.$html_id.'">
                        <div class="card">
                            <legend class="h2 card-header">'.$fa_icon.$param_id.'</legend>
                            <div class="card-body">
                                '.$description.$helptext.$child_parameters.'
                            </div>
                        </div>
                    </fieldset>';
                }
            } else {
                echo build_form_param($param_id, $param, @in_array($param_id, $cache['schema']['required']));
            }
        }
        ?>
        <div class="mt-5 text-center">
            <button type="submit" class="btn btn-lg btn-primary  btn-launch" data-target="#schema-finished">
                <i class="fad fa-rocket-launch"></i> Launch workflow
            </button>
            <p class="small text-danger mt-2 validation-warning" style="display: none;">Please fix validation errors before launching.</p>
        </div>
    </form>

    <div class="toast" role="alert" aria-live="assertive" aria-atomic="true" id="form_validation_error_toast">
        <div class="toast-header">
            <strong class="mr-auto text-danger">Validation error</strong>
            <button type="button" class="ml-2 mb-1 close" data-dismiss="toast" aria-label="Close">
                <span aria-hidden="true">&times;</span>
            </button>
        </div>
        <div class="toast-body">
            <p>There was a problem validating some of your parameters:</p>
            <ul id="validation_fail_list"></ul>
        </div>
    </div>

<?php } // if $cache

include('../includes/footer.php');
