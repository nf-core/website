<?php
$title = 'Tools';
$subtitle = 'The nf-core companion tool, to help with common tasks.';

$markdown_fn = dirname(__FILE__)."/../../markdown/tools/README.md";
$tools_version_fn = dirname(__FILE__)."/../../tools_version.txt";
$tools_all_versions_fn = dirname(__FILE__)."/../../tools_all_versions.txt";
$md_github_url = ' https://github.com/nf-core/tools/blob/master/README.md';
$md_github_raw_url = 'https://raw.githubusercontent.com/nf-core/tools/dev/README.md';

// Get auth secrets
$config = parse_ini_file(dirname(__FILE__)."/../../config.ini");
define('GH_AUTH', base64_encode($config['github_username'].':'.$config['github_access_token']));
// HTTP header to use on GitHub API GET requests
define('GH_API_OPTS',
  stream_context_create([
    'http' => [
      'method' => 'GET',
      'header' => [
        'User-Agent: PHP',
        'Accept:application/vnd.github.mercy-preview+json', // Needed to get topics (keywords) for now
        'Accept:application/vnd.github.luke-cage-preview+json', // Needed to get protected branch required reviews
        "Authorization: Basic ".GH_AUTH
      ]
    ]
  ])
);

// Fetch readme and releases if readme is not found or more than 24 hours old
if(!file_exists($markdown_fn) || filemtime($markdown_fn) < (time() - 60 * 60 * 24)) {

    // Download the readme and cache
    // Build directories if needed
    if (!is_dir(dirname($markdown_fn))) {
        mkdir(dirname($markdown_fn), 0777, true);
    }
    $md_contents = file_get_contents($md_github_raw_url);
    if($md_contents){
        file_put_contents($markdown_fn, $md_contents);
    }

    // Fetch releases
    $tools_versions_url = 'https://api.github.com/repos/nf-core/tools/releases';
    $tools_versions_json = json_decode(file_get_contents($tools_versions_url, false, GH_API_OPTS), true);
    if($tools_versions_json){
        $versions = [];
        foreach($tools_versions_json as $release){
            $versions[] = $release['tag_name'];
        }
        if(count($versions) > 0){
            file_put_contents($tools_all_versions_fn, implode("\n", array_reverse($versions)));
            file_put_contents($tools_version_fn, $versions[0]);
        }
    }
}

$md_trim_before = '## Table of contents';
include('../../includes/header.php');
include('../../includes/footer.php');
?>
