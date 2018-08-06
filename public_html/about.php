<?php
error_reporting(E_ALL);
ini_set('display_errors', TRUE);
ini_set('display_startup_errors', TRUE);

$title = 'About nf-core';
$subtitle = 'Details about the nf-core project - who is involved and how it was started.';
$markdown_fn = '../markdown/about.md';
$no_print_content = true;
include('../includes/header.php');

// Parse YAML contributors file
require_once("../Spyc.php");
$contributors = spyc_load_file('../nf-core-contributors.yaml');
$contributors_html = '<div class="card-deck">';
foreach($contributors['contributors'] as $c){
    // Start card div
    $contributors_html .= '<div class="card contributor"><div class="card-body">';
    // Header, title
    if(array_key_exists('image_fn', $c)){
        $img_path = 'assets/img/contributors-colour/'.$c['image_fn'];
        if(file_exists($img_path))
            $contributors_html .= '<img class="contributor_logo" title="'.$c['full_name'].'" src="'.$img_path.'">';
    }
    if(array_key_exists('full_name', $c)){
        $contributors_html .= '<h5 class="card-title">';
        if(array_key_exists('url', $c))
            $contributors_html .= ' <a href="'.$c['url'].'" target="_blank">';
        $contributors_html .= $c['full_name'];
        if(array_key_exists('url', $c))
            $contributors_html .= ' </a>';
        $contributors_html .= '</h5>';
    }
    if(array_key_exists('affiliation', $c)){
        $contributors_html .= '<h6 class="card-subtitle mb-2 text-muted">';
        if(array_key_exists('affiliation_url', $c))
            $contributors_html .= '<a href="'.$c['affiliation_url'].'" target="_blank">';
        $contributors_html .= $c['affiliation'];
        if(array_key_exists('affiliation_url', $c))
            $contributors_html .= '</a>';
        $contributors_html .= '</h6>';
    }
    // Description
    if(array_key_exists('description', $c))
        $contributors_html .= '<p class="small text-muted">'.$c['description'].'</p> ';
    // Contact person
    $contributors_html .= '<div class="contributor_contact">';
    if(array_key_exists('contact_email', $c)){
        $contributors_html .= '<a href="mailto:'.$c['contact_email'].'" class="badge badge-light" data-toggle="tooltip" title="Primary contact: '.$c['contact_email'].'"><i class="far fa-envelope"></i> ';
        if(array_key_exists('contact', $c))
            $contributors_html .= $c['contact'];
        else
            $contributors_html .= $c['contact_email'];
        $contributors_html .= '</a> ';
    }
    else if(array_key_exists('contact', $c))
            $contributors_html .= $c['contact'].' ';
    if(array_key_exists('contact_github', $c))
        $contributors_html .= '<a href="https://github.com/'.trim($c['contact_github'], '@').'/" target="_blank" class="badge badge-light" data-toggle="tooltip" title="Primary contact: GitHub @'.trim($c['contact_github'], '@').'"><i class="fab fa-github"></i> '.trim($c['contact_github'], '@').'</a> ';
    $contributors_html .= '</div>';
    // Close card div
    $contributors_html .= '</div></div>';
}
$contributors_html .= '</div>';

echo str_replace('<!-- #### CONTRIBUTORS #### -->', $contributors_html, $content);

include('../includes/footer.php');
?>
