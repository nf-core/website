<?php

$title = '<a href="/'.$pipeline->name.'">nf-core/<br class="d-sm-none">'.$pipeline->name.'</a>';
$subtitle = $pipeline->description;

# Details for parsing markdown file, fetched from Github
# Build the remote file path
if(count($path_parts) > 1){
    if(substr($_SERVER['REQUEST_URI'], -3) == '.md'){
        # Clean up URL by removing .md
        header('Location: '.substr($_SERVER['REQUEST_URI'], 0, -3));
        exit;
    }
    $filename = implode('/', array_slice($path_parts, 1)).'.md';
} else {
    $filename = 'README.md';
    $md_trim_before = '# Introduction';
}
# Build the local and remote file paths based on whether we have a release or not
if($pipeline->last_release !== 0){
  $git_branch = 'master';
  $local_fn_base = dirname(dirname(__FILE__))."/markdown/pipelines/".$pipeline->name."/".$pipeline->last_release."/";
} else {
  $git_branch = 'dev';
  $local_fn_base = dirname(dirname(__FILE__))."/markdown/pipelines/".$pipeline->name."/".$pipeline->pushed_at."/";
}
$local_md_fn = $local_fn_base.$filename;
$markdown_fn = 'https://raw.githubusercontent.com/'.$pipeline->full_name.'/'.$git_branch.'/'.$filename;

# Check if we have a local copy of the markdown file and fetch if not
if(file_exists($local_md_fn)){
  $markdown_fn = $local_md_fn;
} else {
  # Build directories if needed
  if (!is_dir(dirname($local_md_fn))) {
    mkdir(dirname($local_md_fn), 0777, true);
  }
  $md_contents = file_get_contents($markdown_fn);
  if($md_contents){
    file_put_contents($local_md_fn, $md_contents);
    $markdown_fn = $local_md_fn;
  } else {
    # Edge case: No releases, but dev branch doesn't exist - use master instead
    $git_branch = 'master';
    $markdown_fn = 'https://raw.githubusercontent.com/'.$pipeline->full_name.'/'.$git_branch.'/'.$filename;
    $md_contents = file_get_contents($markdown_fn);
    if($md_contents){
      file_put_contents($local_md_fn, $md_contents);
      $markdown_fn = $local_md_fn;
    }
  }
}

# Get a navigation tree of all markdown files
// Fetch repo commit tree
$api_opts = stream_context_create([
    'http' => [
        'method' => 'GET',
        'header' => [
            'User-Agent: PHP',
            'Accept:application/vnd.github.mercy-preview+json' // Needed to get topics (keywords) for now
        ]
    ]
]);
$gh_tree_url = "https://api.github.com/repos/{$pipeline->full_name}/git/trees/{$git_branch}?recursive=1";
$local_tree_fn = $local_fn_base.'md_tree.json';
# Check if we have a local copy of the markdown file and fetch if not
if(file_exists($local_tree_fn)){
  $gh_tree_json = file_get_contents($local_tree_fn);
} else {
  # Build directories if needed
  if (!is_dir(dirname($local_tree_fn))) {
    mkdir(dirname($local_tree_fn), 0777, true);
  }
  $gh_tree_json = file_get_contents($gh_tree_url, false, $api_opts);
  if(!in_array("HTTP/1.1 200 OK", $http_response_header)){
      var_dump($http_response_header);
      echo("<!-- Could not fetch nf-core repo tree info! $gh_tree_url -->");
  }
  if($gh_tree_json){
    file_put_contents($local_tree_fn, $gh_tree_json);
  }
}
$gh_tree = json_decode($gh_tree_json);
$md_toc_files = [];
foreach($gh_tree->tree as $tfile){
  if(fnmatch('*.md', $tfile->path) &&
      !fnmatch('.github/*', $tfile->path) &&
      !fnmatch('CODE_OF_CONDUCT.md', $tfile->path) &&
      !fnmatch('docs/README.md', $tfile->path) &&
      !fnmatch('CHANGELOG.md', $tfile->path)
    ){
    $md_toc_files[] = $tfile;
  }
}

# Build ToC of available markdown files
// https://stackoverflow.com/a/19141352/713980
function build_folder_structure(&$md_toc, $path_array) {
    if (count($path_array) > 1) {
        if (!isset($md_toc[$path_array[0]])) {
            $md_toc[$path_array[0]] = array();
        }
        build_folder_structure($md_toc[$path_array[0]], array_splice($path_array, 1));
    } else {
        $md_toc[] = $path_array[0];
    }
}
$md_toc = [];
foreach($md_toc_files as $md_file){
  $path_array = explode('/', $md_file->path);
  build_folder_structure($md_toc, $path_array);
}

# Make a HTML list Table of contents
// https://stackoverflow.com/a/10152223/713980
function array2ul($array, $parents = array()) {
    global $pipeline;
    $out = '<ul class="fa-ul">';
    foreach($array as $key => $elem){
        if(!is_array($elem)){
            $elem = str_replace('.md', '', $elem);
            $elem_url = str_replace('//', '/', '/'.$pipeline->name.'/'.implode('/', $parents).'/'.$elem);
            if($elem_url == '/'.$pipeline->name.'/README'){
                $elem_url = '/'.$pipeline->name;
            }
            $class = '';
            if($_SERVER['REQUEST_URI'] == $elem_url){
              $class = 'active';
            }
            $elem_text = ucfirst(str_replace('_', ' ', $elem));
            $out .= '<li><span class="fa-li"><i class="fas fa-minus" style="font-size: 0.2em; vertical-align: 50%;"></i></span><a href="'.$elem_url.'" class="'.$class.'">'.$elem_text.'</a></li>';
        } else {
            $new_parents = $parents;
            $new_parents[] = $key;
            $elem_text = ucfirst(str_replace('_', ' ', $key));
            $out .= '<li><span class="fa-li"><i class="fas fa-caret-right"></i></span>'.$elem_text.array2ul($elem, $new_parents)."</li>";
        }
    }
    $out .= "</ul>";
    return $out;
}
$md_toc_html = array2ul($md_toc);

# Configs to make relative URLs work
$src_url_prepend = 'https://raw.githubusercontent.com/'.$pipeline->full_name.'/'.$git_branch.'/'.implode('/', array_slice($path_parts, 1, -1)).'/';
$href_url_prepend = $pipeline->name.'/'.implode('/', array_slice($path_parts, 1, -1)).'/';

# Styling
$md_content_replace = array(
    array('# nf-core/'.$pipeline->name.': '),
    array('# ')
);
$html_content_replace = array(
    array('<table>'),
    array('<table class="table">')
);

# Header - keywords
$header_html = '<p class="mb-0">';
foreach($pipeline->topics as $keyword){
  $header_html .= '<a href="/pipelines?q='.$keyword.'" class="badge pipeline-topic">'.$keyword.'</a> ';
}
$header_html .= '</p>';

# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
include('../includes/header.php');

# Pipeline subheader

# Get details for the button to the latest release
function rsort_releases($a, $b){
    $t1 = strtotime($a->published_at);
    $t2 = strtotime($b->published_at);
    return $t2 - $t1;
}
$dev_warning = '';
$archived_warning = '';
if(count($pipeline->releases) > 0){
    usort($pipeline->releases, 'rsort_releases');
    $download_bn = '<a href="'.$pipeline->releases[0]->html_url.'" class="btn btn-success btn-lg">Get version '.$pipeline->releases[0]->tag_name.'</a>';
} else {
    $download_bn = '<a href="'.$pipeline->html_url.'" class="btn btn-success btn-lg">See the development code</a>';
    $dev_warning = '<div class="alert alert-danger">This pipeline is currently in development and does not yet have any stable releases.</div>';
}
if($pipeline->archived){
  $archived_warning = '<div class="alert alert-warning">This pipeline has been archived and is no longer being actively maintained.</div>';
}

# Extra HTML for the header - tags and GitHub URL
?>

<div class="mainpage-subheader-heading">
  <div class="container text-center">
    <?php echo $dev_warning.$archived_warning; ?>
    <p><?php echo $download_bn; ?></p>
    <p class="mb-0"><a href="<?php echo $pipeline->html_url; ?>" class="text-dark"><i class="fab fa-github"></i> <?php echo $pipeline->html_url; ?></a></p>
  </div>
</div>
<div class="triangle subheader-triangle-down"></div>

<div class="container main-content">

<?php

# echo '<pre>'.print_r($pipeline, true).'</pre>';
# echo '<pre>'.print_r($gh_tree_json, true).'</pre>';

echo '<div class="pipeline-page-toc">'.$md_toc_html.'</div>';

# Print rendered markdown made in header.php
echo '<div class="rendered-markdown">';
echo $content;
echo '</div>';

include('../includes/footer.php');
