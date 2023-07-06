<?php
//
// update_fontawesome_icons.json
// ---------------------------
// The pipeline_schema_builder.php page has an icon picker for FontAwesome icons
// so that pipeline authors can icon-ise their pipeline parameters.
// We load Font Awesome on the site through a JavaScript FA kit that always
// gives us the latest version of Font Awesome icons.
// This script runs on a cron job to update a JSON file describing the available
// FA icons so that the icon-picker has the full selection.

# parse the JSON file from Font Awesome
$json_url = 'https://raw.githubusercontent.com/FortAwesome/Font-Awesome/master/metadata/icons.json';
$raw_json = file_get_contents($json_url);
$fa_icons = json_decode($raw_json, true);

# create the icons json object, to be put into iconpicker.js
$icons = [];
foreach ($fa_icons as $d => $icon) {
    foreach ($icon['styles'] as $s) {
        # create one object for each style (e.g. "brand", "solid") of an icon
        $icon_class = 'fa' . $s[0] . ' fa-' . $d; # create 'fab' class name for brand icons, "fas" for solid, etc.
        $search_terms = implode(' ', $icon['search']['terms']);
        $icons[$icon_class] = $search_terms;
    }
}

# Write to a JSON file for the website
$results_fn = dirname(__FILE__) . '/public_html/assets/js/fa-icons.json';
$results_json = json_encode($icons, JSON_PRETTY_PRINT) . "\n";
file_put_contents($results_fn, $results_json);
