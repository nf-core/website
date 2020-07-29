<?php
// Build the HTML docs from a pipeline JSON schema.
// Imported by public_html/pipeline.php

require_once('../includes/parse_md.php');

###############
# RENDER JSON SCHEMA PARAMETERS DOCS
###############
if(file_exists($gh_pipeline_schema_fn)){

  $schema = json_decode(file_get_contents($gh_pipeline_schema_fn), TRUE);
  $hidden_params = [];

  function print_param($is_group, $param_id, $param, $is_required=false){
    global $hidden_params;
    // Is this parameter hidden?
    $is_hidden = false;
    if(array_key_exists("hidden", $param) && $param["hidden"]){
      $is_hidden = true;
    }
    if($is_group){
      $is_hidden = true;
      foreach($param['properties'] as $child_param_id => $child_param){
        if(!array_key_exists("hidden", $child_param) || !$child_param["hidden"]){
          $is_hidden = false;
        }
      }
    } else if($is_hidden) {
        $hidden_params[] = $param_id;
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
      $description = parse_md($param['description'])['content'];
    }

    # Help text
    $help_text_btn = '';
    $help_text = '';
    if(array_key_exists("help_text", $param) && strlen(trim($param['help_text'])) > 0){
      $help_text_btn = '
        <button class="btn btn-sm btn-outline-info ml-2 mb-1 mt-1" data-toggle="collapse" href="#'.$param_id.'-help" aria-expanded="false">
          <i class="fas fa-question-circle"></i> Help
        </button>';
      $help_text = '
        <div class="collapse bg-lightgray col-12 schema-docs-help-text p-2 mb-2 small" id="'.$param_id.'-help">
          '.parse_md($param['help_text'])['content'].'</div>';
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
    $row_class = 'align-items-center';
    $id_cols = 'col-10 col-md-5 col-lg-4 col-xl-3 small-h';
    $h_level = 'h3';
    if($is_group){
      $row_class = 'align-items-baseline mt-5 param-docs-row-group';
      $id_cols = 'col pl-0 h2';
      $h_level = 'h2';
    }
    if($is_hidden){
      $row_class .= ' param-docs-hidden collapse';
    }

    # Build row
    return '
    <div class="row param-docs-row border-bottom '.$row_class.'">
      <div class="'.$id_cols.' param-docs-row-id-col">'.add_ids_to_headers('<'.$h_level.'>'.$fa_icon.$h_text.'</'.$h_level.'>', $is_hidden).'</div>
      <div class="col">'.$param_body.'</div>
      <div class="col-auto text-right">'.$help_text_btn.$hidden_btn.implode(' ', $labels).'</div>
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

  // Hidden params
  if(count($hidden_params) > 0){
    $schema_content .= '
    <div class="alert border bg-light small hidden_params_alert collapse show">
      <p>The following uncommon parameters have been hidden: <code>'.implode('</code>, <code>', $hidden_params).'</code></p>
      <p><a href=".param-docs-hidden, .toc .collapse, .hidden_params_alert" data-toggle="collapse" role="button" aria-expanded="false" aria-controls="collapseExample">Click here</a> to show all hidden params.</p>
    </div>';
  }

  $schema_content .= '</div>'; // .schema-docs
}
