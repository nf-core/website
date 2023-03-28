<?php

$docs_md_base = dirname(dirname(__FILE__)) . '/markdown/';

$path_parts = explode('/', $_GET['path']);
# Clear any empty array elements
$path_parts = array_filter($path_parts);

# First - assume this is usage or contributing docs and try to find the source (old docs structure)
$md_fn = $_GET['path'];
if (substr($md_fn, -3) !== '.md') {
    $md_fn .= '.md';
}

# new docs structure
if (strtolower($path_parts[0]) == 'docs') {
    if (file_exists($docs_md_base . $md_fn)) {
        $markdown_fn = $docs_md_base . $md_fn;
        $md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . $md_fn;
        $section = trim($path_parts[1]);
        include 'documentation.php';
        exit();
    }
    # Directory root with an index.md
    if (file_exists($docs_md_base . substr($md_fn, 0, -3) . '/index.md')) {
        $markdown_fn = $docs_md_base . substr($md_fn, 0, -3) . '/index.md';
        $md_github_url =
            'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . substr($md_fn, 0, -3) . '/index.md';
        $href_url_prepend = substr(basename($md_fn), 0, -3) . '/';
        $section = trim($path_parts[1]);
        include 'documentation.php';
        exit();
    }
    # Directory root with a README.md
    if (file_exists($docs_md_base . substr($md_fn, 0, -3) . '/README.md')) {
        $markdown_fn = $docs_md_base . substr($md_fn, 0, -3) . '/README.md';
        $md_github_url =
            'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . substr($md_fn, 0, -3) . '/README.md';
        $href_url_prepend = substr(basename($md_fn), 0, -3) . '/';
        $section = trim($path_parts[1]);
        include 'documentation.php';
        exit();
    }
}

if (file_exists($docs_md_base . $md_fn)) {
    $markdown_fn = $docs_md_base . $md_fn;
    $md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . $md_fn;
    include '../includes/header.php';
    include '../includes/footer.php';
    exit();
}

# is it the new docs location

# is it a module?
if (strtolower($path_parts[0]) == 'modules') {
    require_once '../includes/functions.php';
    $modules = get_modules();
    $module_name = str_replace('.php', '', $path_parts[1]);

    foreach ($modules as $idx => $module) {
        if (strtolower($module['name']) == strtolower($module_name)) {
            # If capitalization is wrong, redirect because I'm fussy
            if ($module['name'] != $module_name) {
                header('Location: /' . str_replace($path_parts[1], $module['name'], $_GET['path']));
            }
            # Include the script that renders the pipeline page, then exit
            include '../includes/module_page/_index.php';
            exit();
        }
    }
}

# Wasn't docs nor modules - is it a pipeline name?
$pipelines_json = json_decode(file_get_contents('pipelines.json'));

foreach ($pipelines_json->remote_workflows as $pipeline) {
    if (strtolower($pipeline->name) == strtolower($path_parts[0])) {
        # If capilitilsation is wrong, redirect because I'm fussy
        if ($pipeline->name != $path_parts[0]) {
            header('Location: /' . str_replace($path_parts[0], $pipeline->name, $_GET['path']));
        }

        # Include the script that renders the pipeline page, then exit
        include '../includes/pipeline_page/_index.php';
        exit();
    }
}

# Got this far - must be a 404
header('HTTP/1.1 404 Not Found');
include '404.php';
die();
