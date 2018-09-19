<?php
error_reporting(E_ALL);
ini_set('display_errors', TRUE);
ini_set('display_startup_errors', TRUE);

$title = 'About nf-core';
$subtitle = 'Details about the nf-core project - who is involved and how it was started.';
$markdown_fn = '../markdown/about.md';
$no_print_content = true;
$locations = [];
include('../includes/header.php');

// Parse YAML contributors file
require_once("../Spyc.php");
$contributors = spyc_load_file('../nf-core-contributors.yaml');
$contributors_html = '<div class="card contributors-map-card"><div class="card-body" id="contributors-map"></div></div>';
$contributors_html .= '<div class="card-deck">';
foreach($contributors['contributors'] as $idx => $c){
    // Start card div
    $contributors_html .= '<div class="card contributor card_deck_card"><div class="card-body">';
    // Header, title
    $img_path = '';
    if(array_key_exists('image_fn', $c)){
        $img_path = 'assets/img/contributors-colour/'.$c['image_fn'];
        if(file_exists($img_path))
            $contributors_html .= '<img class="contributor_logo" title="'.$c['full_name'].'" src="'.$img_path.'">';
        else $img_path = '';
    }
    $card_id = $idx;
    if(array_key_exists('full_name', $c)){
        $card_id = preg_replace('/[^a-z]+/', '-', strtolower($c['full_name']));
        $contributors_html .= '<h5 class="card-title" id="'.$card_id.'">';
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
    if(array_key_exists('twitter', $c)) {
        if (preg_match("/.*?\/\/twitter.com\/(.*)/", $c['twitter'], $match) && count($match) > 0) {
            $twitter_name = $match[1];
            $contributors_html .= '<a href="https://twitter.com/'.$twitter_name.'/" target="_blank" class="badge badge-light" data-toggle="tooltip" title="@'.$twitter_name.'"><i class="fab fa-twitter"></i> '.$twitter_name.'</a> ';
        }
    }
    $contributors_html .= '</div>';
    // Close card div
    $contributors_html .= '</div></div>';

    // Location JSON
    if(array_key_exists('location', $c)){
        $location['location'] = $c['location'];
        $location['full_name'] = array_key_exists('full_name', $c) ? $c['full_name'] : '';
        $location['card_id'] = $card_id;
        if($img_path){
            $location['image'] = '<br><a href="#'.$card_id.'"><img class="contributor_map_logo" title="'.$location['full_name'].'" src="'.$img_path.'"></a>';
        } else $location['image'] = '';
        array_push($locations, $location);
    }
}
$contributors_html .= '</div>';
$contributors_html .= '<script type="text/javascript">var locations = '.json_encode($locations).'</script>';

echo str_replace('<!-- #### CONTRIBUTORS #### -->', $contributors_html, $content);

include('../includes/footer.php');
?>
<script>
    $(document).ready(function(){
        var latitude = 0.0, longitude = 0.0;
        var map = L.map('contributors-map', {
            zoom: 2
        });
        var greenIcon = new L.Icon({
            iconUrl: 'assets/img/marker-icon-2x-green.png',
            shadowUrl: 'assets/img/marker-shadow.png',
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41]
        });

        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }).addTo(map);

        locations.forEach(function(marker) {
            if (marker != null) {
                latitude += marker.location[0];
                longitude += marker.location[1];
                L.marker(marker.location, {icon: greenIcon}).addTo(map).bindPopup('<a href="#'+marker.card_id+'">'+marker.full_name+'</a>'+marker.image);
            }
        });

        var center = [ latitude / locations.length, longitude / locations.length ];
        map.setView(center);
    });
</script>
