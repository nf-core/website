<?php
// Build the HTML docs from a pipeline JSON schema.
// Imported by public_html/pipeline.php

require_once '../includes/parse_md.php';

###############
# RENDER JSON SCHEMA PARAMETERS DOCS
###############
if (file_exists($gh_pipeline_schema_fn)) {
    $schema = json_decode(file_get_contents($gh_pipeline_schema_fn), true);
    $hidden_params = [];

    function print_param($is_group, $param_id, $param, $is_required = false) {
        global $hidden_params;
        // Is this parameter hidden?
        $is_hidden = false;
        if (array_key_exists('hidden', $param) && $param['hidden']) {
            $is_hidden = true;
        }
        if ($is_group) {
            $is_hidden = true;
            foreach ($param['properties'] as $child_param_id => $child_param) {
                if (!array_key_exists('hidden', $child_param) || !$child_param['hidden']) {
                    $is_hidden = false;
                }
            }
        } elseif ($is_hidden) {
            $hidden_params[] = $param_id;
        }

        # Build heading
        if (array_key_exists('fa_icon', $param)) {
            $fa_icon = '<i class="' . $param['fa_icon'] . ' fa-fw"></i> ';
        } else {
            $fa_icon = '<i class="fad fa-circle fa-fw text-muted"></i> ';
        }
        $h_text = $param['title'];
        if (!$is_group) {
            $h_text = '<code>' . $param_id . '</code>';
        }

        # Description
        $description = '';
        if (array_key_exists('description', $param)) {
            $description = parse_md($param['description'])['content'];
        }
        # parameter type
        $type = '';
        if (array_key_exists('type', $param) && strlen(trim($param['type'])) > 0 && $param['type'] != 'object') {
            $type = is_string($param['type']) ? "'" . $param['type'] . "'" : $param['type'];
            $type =
                '<span class="text-small overflow-scroll w-100"><span class="text-muted">type: </span>' .
                $type .
                '</span>';
        }
        # default value
        $default_val = '';
        if (array_key_exists('default', $param) && strlen(trim($param['default'])) > 0) {
            $default_val = is_string($param['default']) ? "'" . $param['default'] . "'" : $param['default'];
            $default_val =
                '<code class="text-small overflow-scroll w-100"><span class="text-muted">default: </span>' .
                $default_val .
                '</code>';
        }

        # pattern value
        $pattern_val = '';
        if (array_key_exists('pattern', $param) && strlen(trim($param['pattern'])) > 0) {
            $pattern_val =
                '<code class="text-small"><span class="text-muted">pattern: </span>' .
                trim($param['pattern']) .
                '</code>';
        }

        # min/max value
        $minmax_val = '';
        if (array_key_exists('min', $param) && !empty($param['min'])) {
            $minmax_val = '<code class="text-small"><span class="text-muted">min: </span>' . $param['min'] . '</code>';
        } elseif (array_key_exists('max', $param) && !empty($param['max'])) {
            $minmax_val .=
                '<code class="text-small"><span class="text-muted"> max: </span>' . $param['max'] . '</code>';
        }

        # enum value
        $enum_val = '';
        if (array_key_exists('enum', $param) && count($param['enum']) > 0) {
            if (count($param['enum']) == 1) {
                $enum_val =
                    '<code class="text-small"><span class="text-muted">Options: </span>' .
                    implode(' ', $param['enum']) .
                    '</code>';
            } else {
                $param_default = is_string($param['default']) ? "'" . $param['default'] . "'" : $param['default'];
                $dropdown_label =
                    $default_val != ''
                        ? 'Options: <code class="text-small">' . $param_default . '</code> (default)'
                        : 'Options';
                $enum_val =
                    '<div class="btn-group">
                      <button type="button" class="btn btn-light border-0 dropdown-toggle" data-bs-toggle="dropdown" aria-haspopup="true" aria-expanded="false">' .
                    $dropdown_label .
                    '</button>
                    <div class="dropdown-menu">';

                foreach ($param['enum'] as $key => $param_enum) {
                    $enum_val .= '<a class="dropdown-item"><code class="text-small">';
                    $enum_val .= is_string($param_enum) ? "'" . $param_enum . "'" : $param_enum;
                    if ($param_enum == $param['default']) {
                        $enum_val .= ' (default)';
                        $default_val = '';
                    }
                    $enum_val .= '</code></a>';
                }
                $enum_val .= '</div></div>';
            }
        }

        # Help text
        $help_text_btn = '';
        $help_text = '';
        if (array_key_exists('help_text', $param) && strlen(trim($param['help_text'])) > 0) {
            $help_text_btn =
                '
        <button class="btn btn-sm btn-outline-info ms-2 mb-1 mt-1" data-bs-toggle="collapse" href="#' .
                $param_id .
                '-help" aria-expanded="false">
          <i class="fas fa-question-circle"></i> Help
        </button>';
            $help_text =
                '
        <div class="collapse bg-lightgray col-12 param-docs-help-text p-2 mb-2 small d-print-block" id="' .
                $param_id .
                '-help">
          ' .
                parse_md($param['help_text'])['content'] .
                $pattern_val .
                $minmax_val .
                '</div>';
        }
        # Labels
        $labels = [];
        if ($is_required) {
            $labels[] = '<span class="badge bg-warning text-dark small ms-2">required</span>';
        }
        if ($is_hidden) {
            $labels[] = '<span class="badge bg-light text-dark small ms-2">hidden</span>';
        }
        $labels_helpbtn = '';
        if (count($labels) > 0 || strlen($help_text_btn)) {
            $labels_helpbtn =
                '<div class="param_labels_helpbtn d-flex align-items-center">' .
                $help_text_btn .
                implode(' ', $labels) .
                '</div>';
        }

        # Extra group classes
        $row_class = 'align-items-center';
        $id_cols = 'col-10 col-md-5 col-lg-4 col-xl-3 small-h ps-3';
        $h_level = 'h3';
        $description_class = 'small';
        if ($is_group) {
            $row_class = 'align-items-baseline mt-5 param-docs-row-group';
            $id_cols = 'col-12';
            $h_level = 'h2';
            $description_class = 'lead mb-3';
        }
        if ($is_hidden) {
            $row_class .= ' param-docs-hidden collapse d-print-block';
        }

        # Build row
        return '
    <div class="row param-docs-row border-bottom ' .
            $row_class .
            '">
      <div class="' .
            $id_cols .
            ' param-docs-row-id-col align-items-start">' .
            add_ids_to_headers('<' . $h_level . '>' . $fa_icon . $h_text . '</' . $h_level . '>', $is_hidden) .
            $type .
            '</div>
      <div class="col me-auto">' .
            $param_body .
            '</div>
      <div class="col-auto d-flex flex-column align-items-end my-1">' .
            $default_val .
            $enum_val .
            $labels_helpbtn .
            '</div>' .
            $help_text .
            '
    </div>';
    }

    $schema_content = '<div class="param-docs">';
    $schema_content .= _h1('Parameters');
    // Definition groups
    if (isset($schema['allOf']) && count($schema['allOf']) > 0) {
        foreach ($schema['allOf'] as $allof) {
            if (!isset($allof['$ref']) || !isset($schema['definitions'])) {
                continue;
            }
            $group_id = substr($allof['$ref'], 14);
            if (!isset($schema['definitions'][$group_id]) || count($schema['definitions'][$group_id]) == 0) {
                continue;
            }
            $group = $schema['definitions'][$group_id];
            $schema_content .= print_param(true, $group_id, $group);
            // Group-level params
            foreach ($group['properties'] as $param_id => $param) {
                $is_required = false;
                if (array_key_exists('required', $group) && in_array($param_id, $group['required'])) {
                    $is_required = true;
                }
                $schema_content .= print_param(false, '--' . $param_id, $param, $is_required);
            }
        }
    }
    // Top-level params
    if (isset($schema['properties']) && count($schema['properties']) > 0) {
        foreach ($schema['properties'] as $param_id => $param) {
            $is_required = false;
            if (array_key_exists('required', $schema) && in_array($param_id, $schema['required'])) {
                $is_required = true;
            }
            $schema_content .= print_param(false, '--' . $param_id, $param, $is_required);
        }
    }

    // Hidden params
    if (count($hidden_params) > 0) {
        $schema_content .=
            '
    <div class="alert border bg-light small hidden_params_alert collapse show">
      <p>The following uncommon parameters have been hidden: <code>' .
            implode('</code>, <code>', $hidden_params) .
            '</code></p>
      <p><a href=".param-docs-hidden, .toc .collapse, .hidden_params_alert" data-bs-toggle="collapse" role="button" aria-expanded="false" aria-controls="collapseExample">Click here</a> to show all hidden params.</p>
    </div>';
    }

    $schema_content .= '</div>'; // .param-docs
}
