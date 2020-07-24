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

    $div_start = '<div class="param-docs param-docs-'.$level.$hidden_class.' clearfix">';
    $div_end = '</div>';

    # Build heading
    if(array_key_exists("fa_icon", $param)){
      $fa_icon = '<i class="'.$param['fa_icon'].' fa-fw"></i> ';
    } else {
      $fa_icon = '<i class="fad fa-circle fa-fw text-muted"></i> ';
    }
    $h_text = $param_id;
    if($level > 2){ $h_text = '<code>'.$h_text.'</code>'; }
    $heading = _h($level, $fa_icon.$h_text, $is_hidden);

    # Description
    $description = '';
    if(array_key_exists("description", $param)){
      $description = parse_md($param['description']);
    }

    # Help text
    $help_text_btn = '';
    $help_text = '';
    if(array_key_exists("help_text", $param) && strlen(trim($param['help_text'])) > 0){
      // Show help if hidden, as we toggle the whole content
      if(!$is_hidden){
        $help_text_btn = '
          <button class="btn btn-sm btn-outline-info float-right ml-2" data-toggle="collapse" href="#'.$param_id.'-help" aria-expanded="false">
            <i class="fas fa-question-circle"></i> Help
          </button>';
      }
      $help_collapse = $is_hidden ? '' : 'collapse';
      $help_text = '
        <div class="'.$help_collapse.' card schema-docs-help-text" id="'.$param_id.'-help">
          <div class="card-body small text-muted">'.parse_md($param['help_text']).'</div>
        </div>';
    }

    # Hidden button
    $hidden_btn = '';
    if($is_hidden){
      $hidden_btn = '
        <button class="btn btn-sm btn-outline-secondary float-right ml-2" data-toggle="collapse" href="#'.$param_id.'-body" aria-expanded="false">
          <i class="fas fa-eye-slash"></i> Show hidden
        </button>';
    }

    # Labels
    $labels = [];
    if($is_required){
      $labels[] = '<span class="badge badge-warning ml-2">required</span>';
    }
    $labels_helpbtn = '';
    if(count($labels) > 0 || strlen($help_text_btn) || strlen($hidden_btn)){
      $labels_helpbtn = '<div class="param_labels_helpbtn">'.$help_text_btn.$hidden_btn.implode(' ', $labels).'</div>';
    }

    # Body
    $body_collapse = $is_hidden ? 'param-docs-body-hidden collapse' : '';
    $param_body = '<div id="'.$param_id.'-body" class="param-docs-body '.$body_collapse.'">'.$description.$help_text.'</div>';

    return $div_start.$heading.$labels_helpbtn.$param_body.$div_end;
  }

  $schema_content = '<div class="schema-docs">';
  $schema_content .= _h1('Parameters');
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
        $schema_content .= print_param(3, '--'.$child_param_id, $child_param, $is_required);
      }
    }
    // Top-level params
    else {
      $is_required = false;
      if(array_key_exists("required", $schema) && in_array($param_id, $schema["required"])){
        $is_required = true;
      }
      $schema_content .= print_param(3, '--'.$param_id, $param, $is_required);
    }
  }
  $schema_content .= '</div>'; // .schema-docs
}
