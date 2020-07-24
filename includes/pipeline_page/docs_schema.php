<?php
// Build the HTML docs from a pipeline JSON schema.
// Imported by public_html/pipeline.php

require_once('../includes/parse_md.php');

###############
# RENDER JSON SCHEMA PARAMETERS DOCS
###############
if(file_exists($gh_pipeline_schema_fn)){

  $schema = json_decode(file_get_contents($gh_pipeline_schema_fn), TRUE);

  function print_param($level, $param_id, $param, $is_required=false){
    $is_hidden = false;
    if(array_key_exists("hidden", $param) && $param["hidden"]){
      $is_hidden = true;
    }

    # Build heading
    if(array_key_exists("fa_icon", $param)){
      $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw"></i> ';
    } else {
      $fa_icon = '<i class="fad fa-circle fa-fw text-muted"></i> ';
    }
    $h_text = $param_id;
    if($level > 2){ $h_text = '<code>'.$h_text.'</code>'; }
    $heading = _h($level, $fa_icon.$h_text);

    # Description
    $description = '';
    if(array_key_exists("description", $param)){
      $description = parse_md($param['description']);
    }

    # Help text
    $help_text_btn = '';
    $help_text = '';
    if(array_key_exists("help_text", $param)){
      $help_text_btn = '
        <button class="btn btn-sm btn-outline-secondary float-right ml-2" data-toggle="collapse" href="#'.$param_id.'-help" aria-expanded="false">
          <i class="fas fa-question-circle"></i> Help
        </button>';
      $help_text = '
        <div class="collapse card schema-docs-help-text" id="'.$param_id.'-help">
          <div class="card-body small text-muted">'.parse_md($param['help_text']).'</div>
        </div>';
    }

    # Labels
    $labels = [];
    if($is_hidden){
      $labels[] = '<span class="badge badge-secondary ml-2">hidden</span>';
    }
    if($is_required){
      $labels[] = '<span class="badge badge-warning ml-2">required</span>';
    }
    if(count($labels) > 0){
      $labels_str = '<div class="float-right">'.implode(' ', $labels).'</div>';
    }

    return $labels_str.$heading.$help_text_btn.$description.$help_text;
  }

  $schema_content = '<div class="schema-docs">';
  $schema_content .= _h(1, 'Parameters');
  foreach($schema["properties"] as $param_id => $param){
    // Groups
    if($param["type"] == "object"){
      $schema_content .= print_param(2, $param_id, $param);
      // Group-level params
      foreach($param["properties"] as $child_param_id => $child_param){
        $is_required = false;
        if(array_key_exists("required", $param) && in_array($child_param_id, $param["required"])){
          $is_required = true;
        }
        $schema_content .= print_param(3, $child_param_id, $child_param, $is_required);
      }
    }
    // Top-level params
    else {
      $is_required = false;
      if(array_key_exists("required", $schema) && in_array($param_id, $schema["required"])){
        $is_required = true;
      }
      $schema_content .= print_param(3, $param_id, $param, $is_required);
    }
  }
  $schema_content .= '</div>'; // .schema-docs
}
