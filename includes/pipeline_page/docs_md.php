<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

# Details for parsing markdown file, fetched from Github
# Build the remote file path
# Special case - root docs is allow
# General docs page

# Build the local and remote file paths based on whether we have a release or not
if ($release !== 'dev') {
    $local_fn_base =
        dirname(dirname(dirname(__FILE__))) . '/markdown/pipelines/' . $pipeline->name . '/' . $release . '/';
} else {
    $local_fn_base =
        dirname(dirname(dirname(__FILE__))) .
        '/markdown/pipelines/' .
        $pipeline->name .
        '/' .
        $pipeline->pushed_at .
        '/';
}
$local_md_fn = $local_fn_base . $filename;
$markdown_fn = 'https://raw.githubusercontent.com/' . $pipeline->full_name . '/' . $release . '/' . $filename;

# Check if we have a local copy of the markdown file and fetch if not
if (file_exists($local_md_fn)) {
    $markdown_fn = $local_md_fn;
} else {
    # Build directories if needed
    if (!is_dir(dirname($local_md_fn))) {
        mkdir(dirname($local_md_fn), 0777, true);
    }
    $md_contents = file_get_contents($markdown_fn);
    if ($md_contents) {
        file_put_contents($local_md_fn, $md_contents);
        $markdown_fn = $local_md_fn;
    }
}

# Configs to make relative URLs work
$href_url_prepend =  '/' . $release . '/' . dirname($filename) . '/';
$href_url_prepend = preg_replace('/\/\/+/', '/', $href_url_prepend);
$href_url_prepend = 'https://github.com/' . $pipeline->full_name . '/blob/' . $href_url_prepend;

$src_url_prepend = '/' . $pipeline->name . '/' . $release . '/' . dirname($filename) . '/';
$src_url_prepend= preg_replace('/\/\/+/', '/', $src_url_prepend);
$src_url_prepend = 'https://raw.githubusercontent.com/sanger-tol' . $src_url_prepend;

$href_url_suffix_cleanup = '\.md';

# Markdown cleanup
$md_content_replace[] = ['/# sanger-tol\/' . $pipeline->name . ': /', '# '];
$md_content_replace[] = [
    '/# !\[sanger-tol\/' . $pipeline->name . '\]\(images\/sanger-tol-' . $pipeline->name . '_logo.png\)/',
    '',
];
$md_content_replace[] = ['/(## :warning:)(.*?)( files\._)/s', ''];

// Footer source link
$md_github_url = 'https://github.com/' . $pipeline->full_name . '/blob/' . $release . '/' . $filename;
