<?php
// Find a search term - return results for HTML page or API

require '../vendor/autoload.php';
use Spyc;

// $search_term - should be available from include
$search_results = [
    'pipelines' => [],
    'documentation' => [],
];

//
// Search pipeline names and keywords
//
$pipelines_json = json_decode(file_get_contents('pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
// Sort alphabetically
usort($pipelines, function ($a, $b) {
    strcmp($a->name, $b->name);
});

foreach ($pipelines as $pipeline) {
    $hit_name = false;
    $hit_keywords = false;
    $hit_description = false;
    // Search pipeline name
    $pipeline->url = '/' . $pipeline->name;
    if (stripos($pipeline->name, $search_term) !== false) {
        $hit_name = true;
        $pipeline->full_name = preg_replace("/($search_term)/i", "<mark>$1</mark>", $pipeline->full_name);
        $pipeline->name = preg_replace("/($search_term)/i", "<mark>$1</mark>", $pipeline->name);
    }
    // Search pipeline keywords
    foreach ($pipeline->topics as $idx => $kw) {
        if (stripos($kw, $search_term) !== false) {
            $hit_keywords = true;
            $pipeline->topics[$idx] = preg_replace("/($search_term)/i", "<mark>$1</mark>", $kw);
        }
    }
    // Search pipeline description
    if (stripos($pipeline->description, $search_term) !== false) {
        $hit_description = true;
        $pipeline->description = preg_replace("/($search_term)/i", "<mark>$1</mark>", $pipeline->description);
    }
    if ($hit_name || $hit_keywords || $hit_description) {
        $search_results['pipelines'][] = [
            'pipeline' => $pipeline,
            'hit_name' => $hit_name,
            'hit_keywords' => $hit_keywords,
            'hit_description' => $hit_description,
        ];
    }
}
// Sort by hits - name is best, keyword only second, description only third
usort($search_results['pipelines'], function ($a, $b) {
    if ($a['hit_name'] && $b['hit_name']) {
        return 0;
    }
    if ($a['hit_name']) {
        return -1;
    }
    if ($b['hit_name']) {
        return 1;
    }
    if ($a['hit_keywords'] && $b['hit_keywords']) {
        return 0;
    }
    if ($a['hit_keywords']) {
        return -1;
    }
    if ($b['hit_keywords']) {
        return 1;
    }
    if ($a['hit_description'] && $b['hit_description']) {
        return 0;
    }
    if ($a['hit_description']) {
        return -1;
    }
    if ($b['hit_description']) {
        return 1;
    }
});

//
// Search markdown files
//
$md_base = dirname(dirname(__FILE__)) . '/markdown';
foreach (['usage', 'contributing'] as $docs_type) {
    foreach (glob("$md_base/$docs_type/*.md") as $file) {
        $content = file_get_contents("$file");
        $match_pos = stripos($content, $search_term);
        if ($match_pos !== false) {
            // Highlight the match
            $content = preg_replace("/($search_term)/i", "<mark>$1</mark>", $content);
            // Get an excerpt around the first match
            $match_string =
                '&hellip;' . substr($content, $match_pos - 40 - 6, strlen($search_term) + 80 + 7 + 6) . '&hellip;';

            // Find the title of the docs page
            $title = ucfirst(basename($file, '.md'));
            $subtitle = '';
            $md = $content;
            if (substr($content, 0, 3) == '---') {
                $md_parts = explode('---', $content, 3);
                if (count($md_parts) == 3) {
                    $meta = spyc_load($md_parts[1]);
                    $md = $md_parts[2];
                    if (isset($meta['title'])) {
                        $title = $meta['title'];
                    }
                    if (isset($meta['subtitle'])) {
                        $subtitle = $meta['subtitle'];
                    }
                }
            }

            // Flags for where this was found
            $hit_title = stripos($meta['title'], $search_term);
            $hit_subtitle = stripos($meta['subtitle'], $search_term);
            $hit_content = stripos($md, $search_term);

            $search_results['documentation'][] = [
                'match_string' => $match_string,
                'path' => "$file",
                'url' => "/$docs_type/" . basename($file, '.md'),
                'title' =>
                    ucfirst($docs_type) .
                    ' &raquo; ' .
                    (isset($meta['title']) ? $meta['title'] : ucfirst(basename($file, '.md'))),
                'short_title' => $title,
                'subtitle' => $subtitle,
                'hit_title' => $hit_title,
                'hit_subtitle' => $hit_subtitle,
                'hit_content' => $hit_content,
            ];
        }
    }
}

// Sort by hits - name is best, subtitle second, contents third
usort($search_results['documentation'], function ($a, $b) {
    if ($a['hit_title'] && $b['hit_title']) {
        return 0;
    }
    if ($a['hit_title']) {
        return -1;
    }
    if ($b['hit_title']) {
        return 1;
    }
    if ($a['hit_subtitle'] && $b['hit_subtitle']) {
        return 0;
    }
    if ($a['hit_subtitle']) {
        return -1;
    }
    if ($b['hit_subtitle']) {
        return 1;
    }
    if ($a['hit_content'] && $b['hit_content']) {
        return 0;
    }
    if ($a['hit_content']) {
        return -1;
    }
    if ($b['hit_content']) {
        return 1;
    }
});
