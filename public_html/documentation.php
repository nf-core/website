<?php
require_once '../includes/functions.php';
require_once '../includes/parse_md.php';

########
## Configure page header
########
$content = parse_md($markdown_fn)['content'];

$no_auto_toc = true;

# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
$sidebar_nav_elements = [];

# function to generate an array of the file structure
function dir_tree($dir) {
    global $docs_md_base;
    $files = array_map('basename', glob("$dir/*"));
    foreach ($files as $file) {
        if (is_dir("$dir/$file")) {
            $tree[$file] = dir_tree("$dir/$file");
        } else {
            $fm = file_get_contents("$dir/$file");
            $fm = parse_md_front_matter($fm);
            $tree[] = [
                'title' => $fm['meta']['title'],
                'weight' => $fm['meta']['menu']['main']['weight'] ? $fm['meta']['menu']['main']['weight'] : 100,
                'url' => str_replace($docs_md_base, '', preg_replace('/\.md$/', '', "$dir/$file")),
            ];
        }
        # have items with lower weight come first, rest sort by title
        array_multisort(array_column($tree, 'weight'), SORT_ASC, array_column($tree, 'title'), SORT_ASC, $tree);
    }
    return $tree;
}

$sidebar_nav_elements = dir_tree($docs_md_base . 'docs');
//Remove landing page from navigation
unset($sidebar_nav_elements['0']);
krsort($sidebar_nav_elements, SORT_ASC); # sort Usage before Contributing
# build html for the sidebar nav

function build_sidebar_nav($elements) {
    global $sidebar_nav;
    global $md_fn;
    foreach ($elements as $name => $element) {
        if (!isset($element['url'])) {
            $path = explode('/', $element[0]['url']);
            $show = strpos($md_fn, implode('/', array_slice($path, 0, -1))) !== false ? 'show' : '';
            $is_open = $show == 'show' ? 'true' : 'false';
            $id = str_replace(' ', '-', strtolower(implode('-', array_slice($path, -3, 2))));
            $id = str_replace(':', '', $id);
            $sidebar_nav .=
                '<button class="btn d-inline-flex align-items-center" data-bs-toggle="collapse" data-bs-target="#' .
                $id .
                '" aria-expanded="' .
                $is_open .
                '" aria-current="' .
                $is_open .
                '">
                            <i class="fas fa-angle-right me-3"></i><strong>' .
                ucwords($name) .
                '
                        </strong></button>';
            $sidebar_nav .=
                '<nav class="collapse ' . $show . '" id="' . $id . '"><ul class="list-unstyled fw-normal ps-3 small">';
            build_sidebar_nav($element);
            $sidebar_nav .= '</ul></nav>';
        } else {
            $active = '';
            if (
                $md_fn == $element['url'] . '.md' ||
                substr($md_fn, 0, -3) . '/index' == $element['url'] ||
                substr($md_fn, 0, -3) . '/README' == $element['url']
            ) {
                $active = 'active';
            }
            $sidebar_nav .=
                '<li><a href="/' .
                $element['url'] .
                '"  class="d-inline-flex align-items-center py-1 px-2 ' .
                $active .
                '">' .
                $element['title'] .
                '</a></li>';
        }
    }
}

$sidebar_nav = '<nav class="sidebar-nav side-sub-subnav sticky-top"><ul class="ps-0 d-flex flex-column">';
$sidebar_nav .= build_sidebar_nav($sidebar_nav_elements);
$sidebar_nav .= '</ul></nav>';

# ToC
$toc_nav = '<nav class="toc auto-toc mt-2 flex-column border-start">';
$toc_nav .= generate_toc($content);
$toc_nav .=
    '<p class="small text-end mt-3 d-none d-lg-block"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
$toc_nav .= '</nav>';

$md_content_replace[] = ['<!-- usage_toc -->'];
if (in_array($_GET['path'], ['docs'])) {
    $inline_toc = '<div class="row"><div class="col-12 col-md-6"><div class="mb-3 h3">Usage</div><ul>';
    foreach ($sidebar_nav_elements['usage'] as $mdfile => $mdcontent) {
        if (isset($mdcontent['url'])) {
            $mdcontent['url'] = str_replace('docs/', '', $mdcontent['url']);
            $inline_toc .= '<li>' . '<a href="' . $mdcontent['url'] . '">' . $mdcontent['title'] . '</li>';
        } else {
            $mdcontent[0]['url'] = str_replace('docs/', '', $mdcontent[0]['url']);
            $inline_toc .= '<li>' . ucfirst($mdfile) . '</li>';
            $inline_toc .= '<ul>';
            foreach ($mdcontent as $dropfile => $dropcontent) {
                $dropcontent['url'] = str_replace('docs/', '', $dropcontent['url']);
                $inline_toc .= '<li><a href="' . $dropcontent['url'] . '">' . $dropcontent['title'] . '</li>';
            }
            $inline_toc .= '</ul>';
        }
    }
    $inline_toc .= '</ul></div>';
    $inline_toc .= '<div class="col-12 col-md-6"><div class="mb-3 h3">Contributing</div><ul>';
    foreach ($sidebar_nav_elements['contributing'] as $mdfile => $mdcontent) {
        if (isset($mdcontent['url'])) {
            $mdcontent['url'] = str_replace('docs/', '', $mdcontent['url']);
            $inline_toc .= '<li>' . '<a href="' . $mdcontent['url'] . '">' . $mdcontent['title'] . '</li>';
        } else {
            $mdcontent[0]['url'] = str_replace('docs/', '', $mdcontent[0]['url']);
            $inline_toc .= '<li>' . ucfirst($mdfile) . '</li>';
            $inline_toc .= '<ul>';
            foreach ($mdcontent as $dropfile => $dropcontent) {
                $dropcontent['url'] = str_replace('docs/', '', $dropcontent['url']);
                if (is_array($dropcontent[0])) {
                    $dropcontent[0]['url'] = str_replace('docs/', '', $dropcontent[0]['url']);
                    $inline_toc .= '<ul>';
                    $inline_toc .=
                        '<li>' . '<a href="' . $dropcontent[0]['url'] . '">' . $dropcontent[0]['title'] . '</li>';
                    $inline_toc .= '</ul>';
                } else {
                    $inline_toc .= '<li>' . '<a href="' . $dropcontent['url'] . '">' . $dropcontent['title'] . '</li>';
                }
            }
            $inline_toc .= '</ul>';
        }
    }
    $inline_toc .= '</ul></div></div>';
    $md_content_replace[] = ['/<!-- inline_toc -->/', $inline_toc];
}

include '../includes/header.php';
?>
<div class="container-xxl main-content">
    <?php
    $main_content = '<div class="row">';

    # left sidebar
    $main_content .= '<div class="col-12 col-lg-2">';
    $main_content .= $sidebar_nav;
    $main_content .= '</div>';

    # right sidebar
    $main_content .= '<div class="col-12 col-lg-2 order-lg-last ps-2 h-100 sticky-top"><div class="side-sub-subnav">';
    $main_content .= substr_count($toc_nav, 'nav-link') > 1 ? $toc_nav : '';
    // $main_content .= $toc_nav ;
    $main_content .= '</div></div>'; # end of the sidebar col
    // # main content
    $main_content .= '<div class="col-12 col-lg-8"><div class="rendered-markdown">' . $content . '</div></div>';
    $main_content .= '</div>'; # end of the row
    echo $main_content;
    ?>
</div>

<?php include '../includes/footer.php';
