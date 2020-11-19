<?php

$docs_md_base = dirname(dirname(__file__))."/markdown/";

# First - assume this is usage or developer docs and try to find the source
$md_fn = $_GET['path'];
if(substr($md_fn, -3) !== '.md'){
    $md_fn .= '.md';
}
if(file_exists($docs_md_base.$md_fn)){
    $markdown_fn = $docs_md_base.$md_fn;
    $md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/'.$md_fn;
    include('../includes/header.php');
    include('../includes/footer.php');
    exit;
}


# Wasn't docs - is it a pipeline name?
$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$path_parts = explode('/', $_GET['path']);
# Clear any empty array elements
$path_parts = array_filter($path_parts);
foreach($pipelines_json->remote_workflows as $pipeline){
    if(strtolower($pipeline->name) == strtolower($path_parts[0])){
        # If capilitilsation is wrong, redirect because I'm fussy
        if($pipeline->name != $path_parts[0]){
            header('Location: /'.str_replace($path_parts[0], $pipeline->name, $_GET['path']));
        }

        # Include the script that renders the pipeline page, then exit
        include('../includes/pipeline_page/_index.php');
        exit;
    }
}

# Got this far - must be a 404
header('HTTP/1.1 404 Not Found');
include('404.php');
die();
