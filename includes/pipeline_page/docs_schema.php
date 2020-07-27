<?php
// Build the HTML docs from a pipeline JSON schema.
// Imported by public_html/pipeline.php

require_once('../includes/parse_md.php');

###############
# RENDER JSON SCHEMA PARAMETERS DOCS
###############
if(file_exists($gh_pipeline_schema_fn)){

  $schema = json_decode(file_get_contents($gh_pipeline_schema_fn), TRUE);

  function print_param($is_group, $param_id, $param, $is_required=false){
    $is_hidden = false;
    $hidden_class = '';
    if(array_key_exists("hidden", $param) && $param["hidden"]){
      $is_hidden = true;
      $hidden_class = ' param-docs-hidden';
    }
    if($param['type'] == 'object'){
      $is_hidden = true;
      foreach($param['properties'] as $child_param_id => $child_param){
        if(!array_key_exists("hidden", $child_param) || !$child_param["hidden"]){
          $is_hidden = false;
        }
      }
    }

    # Build heading
    if(array_key_exists("fa_icon", $param)){
      $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw"></i> ';
    } else {
      $fa_icon = '<i class="fad fa-circle fa-fw text-muted"></i> ';
    }
    $h_text = $param_id;
    if(!$is_group){ $h_text = '<code>'.$h_text.'</code>'; }

    # Description
    $description = '';
    if(array_key_exists("description", $param)){
      $description = parse_md($param['description']);
    }

    # Help text
    $help_text_btn = '';
    $help_text = '';
    if(array_key_exists("help_text", $param) && strlen(trim($param['help_text'])) > 0){
      $help_text_btn = '
        <button class="btn btn-sm btn-outline-info float-right ml-2" data-toggle="collapse" href="#'.$param_id.'-help" aria-expanded="false">
          <i class="fas fa-question-circle"></i> Help
        </button>';
      $help_text = '
        <div class="collapse card col-12 schema-docs-help-text" id="'.$param_id.'-help">
          <div class="card-body small text-muted">'.parse_md($param['help_text']).'</div>
        </div>';
    }

    # Labels
    $labels = [];
    if($is_required){
      $labels[] = '<span class="small badge badge-warning ml-2">required</span>';
    }
    if($is_hidden){
      $labels[] = '<span class="small badge badge-light ml-2">hidden</span>';
    }
    $labels_helpbtn = '';
    if(count($labels) > 0 || strlen($help_text_btn)){
      $labels_helpbtn = '<div class="param_labels_helpbtn">'.$help_text_btn.implode(' ', $labels).'</div>';
    }

    # Body
    $param_body = '<div id="'.$param_id.'-body" class="param-docs-body small">'.$description.'</div>';

    # Extra group classes
    $mt = '';
    $id_cols = 'col-10 col-md-5 col-lg-4 col-xl-3 ';
    if($is_group){
      $mt = 'mt-5';
      $id_cols = 'col h2';
    }

    # Build row
    return '
    <div class="row '.$hidden_class.' params-docs-row border-bottom pt-2 '.$mt.'">
      <div class="'.$id_cols.'">'.$fa_icon.$h_text.'</div>
      <div class="col">'.$param_body.'</div>
      <div class="col-auto">'.$help_text_btn.$hidden_btn.implode(' ', $labels).'</div>
      '.$help_text.'
    </div>';
  }

  $schema_content = '<div class="schema-docs">';
  $schema_content .= _h1('Parameters');
  foreach($schema["properties"] as $param_id => $param){
    // Groups
    if($param["type"] == "object"){
      $schema_content .= print_param(true, $param_id, $param);
      // Group-level params
      foreach($param["properties"] as $child_param_id => $child_param){
        $is_required = false;
        if(array_key_exists("required", $param) && in_array($child_param_id, $param["required"])){
          $is_required = true;
        }
        $schema_content .= print_param(false, '--'.$child_param_id, $child_param, $is_required);
      }
    }
    // Top-level params
    else {
      $is_required = false;
      if(array_key_exists("required", $schema) && in_array($param_id, $schema["required"])){
        $is_required = true;
      }
      $schema_content .= print_param(false, '--'.$param_id, $param, $is_required);
    }
  }
  $schema_content .= '</div>'; // .schema-docs
}
