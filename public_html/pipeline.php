<?php

$title = 'nf-core/<br class="d-sm-none">'.$pipeline->name;
$subtitle = $pipeline->description;

# Details for parsing markdown readme, fetched from Github
$markdown_fn = 'https://raw.githubusercontent.com/'.$pipeline->full_name.'/master/README.md';
$md_trim_before = '# Introduction';
$src_url_prepend = 'https://raw.githubusercontent.com/'.$pipeline->full_name.'/master/';
$href_url_prepend = 'https://github.com/'.$pipeline->full_name.'/blob/master/';

# Get details for the button to the latest release
function rsort_releases($a, $b){
    $t1 = strtotime($a->published_at);
    $t2 = strtotime($b->published_at);
    return $t2 - $t1;
}
if(count($pipeline->releases) > 0){
    usort($pipeline->releases, 'rsort_releases');
    $header_btn_url = $pipeline->releases[0]->html_url;
    $header_btn_text = 'Version '.$pipeline->releases[0]->tag_name;
}

# Extra HTML for the header - tags and GitHub URL
$header_html = '
<div class="row">
    <div class="col-sm-6">
        <a href="'.$pipeline->html_url.'" class="text-light"><i class="fab fa-github"></i> '.$pipeline->html_url.'</a>
    </div>
    <div class="col-sm-6 text-right">';
    foreach($pipeline->topics as $keyword){
        $header_html .= '<a href="/pipelines?q='.$keyword.'" class="badge font-weight-light btn btn-sm btn-outline-light">'.$keyword.'</a> ';
    }
$header_html .= '
    </div>
</div>';

include('../includes/header.php');

# echo '<pre>'.print_r($pipeline, true).'</pre>';

include('../includes/footer.php');
