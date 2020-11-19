<?php

$cache_dir = dirname(dirname(__FILE__)).'/api_cache/json_builder';
$post_content_type = 'json_schema';
$post_keys = ['version', 'schema', 'status'];
require_once('../includes/json_schema.php');

// Got this far without printing JSON - build web GUI
$title = 'Parameter schema';
$subtitle = 'Customise a JSON Schema for your pipeline parameters';
if($cache) $import_schema_builder = true;
$mainpage_container = false;
include('../includes/header.php');
?>
<div class="container container-xl">

<button class="btn btn-outline-info float-right ml-2 mb-2" data-toggle="collapse" href="#page_help"><i class="fas fa-question-circle mr-1"></i> Help</button>
<p class="mt-5">nf-core pipelines have a file in their root directories called <code>nextflow_schema.json</code> which
describes the input parameters that the pipeline accepts.
This page helps pipeline authors to build their pipeline schema file by using a graphical interface.</p>

<div class="card collapse mb-2" id="page_help">
  <div class="card-body">
    <h5>nf-core schema</h5>

    <p>nf-core schema files use the <a href="https://json-schema.org/" target="_blank">JSON Schema</a> <em>Draft 7</em> standard:</p>
    <ul>
        <li>Pipeline parameters should be described as <code>properties</code> either in the top-level schema, or in subschema within <code>definitions</code>.</li>
        <li>The <code>definitions</code> subschemas are used to group parameters for the user-interface.
            <ul>
                <li>They are combined in the main schema using <code>allOf</code> for parameter validation.</li>
                <li><code>allOf</code> is a list, the order of this list defines the order that the groups are displayed.</li>
            </ul>
        </li>
        <li>Ungrouped params in the main schema <code>properties</code> are fine, but they will always be sorted at the end, after all <code>definitions</code> groups.</li>
    </ul>
    <p>Although the pipeline parameter validation can handle nesting of parameters in schema <code>objects</code> (eg. <code>params.foo.bar = "baz"</code>) and multiple-levels
        of nesting <code>definitions</code> subschema groups, this is not supported by nf-core - tools such as this builder are likely to behave unpredictably.</p>
    <p>We use a couple of extra JSON keys in addition to the standard JSON Schema set:</p>
    <ul>
        <li><code>help_text</code>, a longer description providing more in-depth help. Typically <code>description</code> is just a few words long and the longer help text is shown when a user requests it.</li>
        <li><code>hidden: True</code>, which tells tools to ignore this <code>param</code> in interfaces by default.</li>
        <li><code>fa_icon</code>, a <a href="https://fontawesome.com/" target="_blank">fontawesome.com</a> icon for use in web interfaces (eg: <code>&lt;i class="fas fa-flask"&gt;&lt;/i&gt;</code> - <i class="fas fa-flask"></i> )</li>
    </ul>

    <h5>Schema tips:</h5>
    <ul>
        <li><code>string</code>, <code>number</code>, <code>integer</code> and <code>range</code> parameters can take a list of <em>enumerated values</em> - a set of allowed values. User interfaces will then display a dropdown select-list.</li>
        <li><code>string</code> params can also have a <em>pattern</em> - a regular expression to validate the input.</li>
        <li><code>range</code> parameters can have either a <em>Minimum</em> or a <em>Maximum</em> or both.</li>
        <li>All of the above settings are accessible through the settings <i class="fas fa-cog"></i></li>
    </ul>

    <h5>Builder tips:</h5>
    <ul>
        <li>Click the <i class="fas fa-cog"></i> icon on the right to access more settings.</li>
        <li>Click and drag the <i class="fas fa-grip-vertical"></i> icon on the left to re-order parameters and groups.</li>
        <li>Click the <i class="fas fa-folder-download"></i> icon to select top-level parameters to move into the group.</li>
        <li>The <i class="fas fa-icons"></i> icon allows you to set a custom icon for the parameter or group.</li>
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
  </div>
</div>

<?php if(!$cache){ ?>

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
<div class="container container-xl main-content">

    <div class="schema-gui-header sticky-top">
        <div class="row align-items-center">
            <div class="col">
                <div class="btn-group mr-3" role="group">
                    <button class="btn btn-outline-secondary add-param-btn"><i class="fas fa-plus-square mr-1"></i> Add parameter</button>
                    <button class="btn btn-outline-secondary add-group-btn"><i class="fas fa-folder-plus mr-1"></i> Add group</button>
                </div>
                <div class="btn-group mr-3" role="group">
                    <button class="btn btn-outline-secondary collapse-groups-btn"><i class="fas fa-folder mr-1"></i> Collapse groups</button>
                    <button class="btn btn-outline-secondary expand-groups-btn"><i class="fas fa-folder-open mr-1"></i> Expand groups</button>
                </div>
                <button class="btn btn-outline-secondary to-top-btn schema-panel-btn" data-target="#schema-builder"><i class="fas fa-arrow-to-top mr-1"></i> Back to top</button>
            </div>
            <div class="col-auto">
                <button class="btn btn-primary schema-panel-btn" data-target="#schema-finished"><i class="fas fa-check-square mr-1"></i> Finished</button>
            </div>
        </div>
    </div>
    <div class="schema-panel" id="schema-builder"></div>
    <div class="schema-panel" id="schema-finished" style="display:none;">
        <div class="card mb-5">
            <div class="card-body">
                <h2 id="schema-send-status">Ok, that's it - done!</h2>
                <p>If you still have <code>nf-core schema build</code> running in the background,
                    it should now update with your new schema. You can close this window.</p>
                <p>If you didn't use <code>nf-core schema build</code> or it has already exited,
                    copy the schema below and paste it in to your pipeline's <code>nextflow_schema.json</code> file.</p>
                <p>Remember that you can run <code>nf-core schema build</code> as many times as you like to make incremental updates.</p>
                <p class="text-center"><button class="btn btn-outline-secondary copy-schema-btn">
                    <i class="far fa-copy mr-2"></i> Copy pipeline schema
                </button><p>
                <p class="text-center"><button class="btn btn-outline-secondary back-to-editor-btn schema-panel-btn" data-target="#schema-builder">
                    <i class="fas fa-brackets-curly mr-2"></i> Back to editor
                </button></p>
            </div>
        </div>
    </div>

    <h3>Pipeline JSON Schema</h3>
    <p>This is the schema for your pipeline. As you change values in the form above, it will update. When you are finished, click <em>Finished</em> in the top toolbar.</p>
    <textarea id="json_schema" class="form-control text-monospace disabled" disabled rows="30"><?php echo json_encode($cache['schema'], JSON_UNESCAPED_SLASHES | JSON_PRETTY_PRINT); ?></textarea>

    <!-- Params schema settings modal -->
    <div class="modal fade" id="settings_modal" tabindex="-1" role="dialog" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <span class="modal-title h4"></span>
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
        <div class="modal-dialog modal-xl" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <span class="modal-title h4"></span>
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
                                        <i class="fab fa-markdown"></i> <a href="https://www.markdownguide.org/cheat-sheet/" target="_blank">Markdown</a> is supported,
                                        but this will be shown raw on the command line so please keep it simple.
                                    </small>
                                </div>
                            </div>
                            <div class="tab-pane fade" id="tab-helptext-preview" role="tabpanel">
                                <p>Command-line:</p>
                                <pre><span class="helptext-cli-preview-title font-weight-bold"></span>
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

    <!-- Moving multiple params into group modal -->
    <div class="modal fade" id="multi_select_modal" tabindex="-1" role="dialog" aria-hidden="true">
        <div class="modal-dialog modal-dialog-scrollable modal-lg" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <h4>Move parameters into "<span></span>"</h4>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                      <span aria-hidden="true">&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <div class="alert alert-info mb-0" id="no_params_alert">No ungrouped parameters available.</div>
                    <div class="params_table">
                        <div class="input-group">
                            <div class="input-group-prepend">
                                <div class="input-group-text" id="btnGroupAddon"><i class="fad fa-search"></i></div>
                            </div>
                            <input type="text" class="form-control" placeholder="Search parameters" id = "search_parameters" aria-label="Search parameters">
                        </div>
                        <p class="small text-muted mb-0 mt-2">Select a range of parameters by holding the shift key while selecting two parameters.</p>
                        <table class="table table-sm" id="params_table">
                            <thead>
                                <tr>
                                <th scope="col" style="width:1%"></th>
                                <th scope="col">Parameter</th>
                                <th scope="col">Description</th>
                                </tr>
                            </thead>
                            <tbody>
                            </tbody>
                        </table>
                    </div>
                </div>
                <div class="modal-footer">
                    <div class="btn-group" role="group">
                        <button type="button" class="btn btn-outline-secondary" id="select_all_params">Select all</button>
                        <button type="button" class="btn btn-outline-secondary" id="deselect_all_params">Deselect all</button>
                    </div>
                    <div class="col text-right pr-0">
                    <button type="button" class="btn btn-secondary" data-dismiss="modal">Cancel</button>
                    <button type="button" class="btn btn-primary disabled" id="move_params">
                        Move parameters
                    </button>
                    </div>
                </div>
            </div>
        </div>
    </div>
</div> <!-- .container-xl -->

<?php } // if $cache

include('../includes/footer.php');
