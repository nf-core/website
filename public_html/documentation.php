<?php
require_once('../includes/functions.php');
require_once('../includes/parse_md.php');

########
## Configure page header
########
$content = parse_md($markdown_fn)['content'];

$no_auto_toc = true;

# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
$files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($docs_md_base . 'docs/', RecursiveDirectoryIterator::FOLLOW_SYMLINKS | RecursiveDirectoryIterator::SKIP_DOTS), RecursiveIteratorIterator::LEAVES_ONLY);
$sidebar_nav_elements = [];
foreach ($files as $name => $object) {
    $md_full = file_get_contents($name);
    if ($md_full !== false) {
        $fm = parse_md_front_matter($md_full);
        $parent = $fm['meta']['menu']['main']['parent'] ? $fm['meta']['menu']['main']['parent'] : end(explode('/', dirname($name)));
        $sidebar_nav_elements[$parent] = [
            'title' => $fm['meta']['title'],
            'parent' => $fm['meta']['menu']['main']['parent'] ? $fm['meta']['menu']['main']['parent'] : end(explode('/',dirname($name))),
            'weight' => $fm['meta']['menu']['main']['weight'] ? $fm['meta']['menu']['main']['weight'] : 0,
            'url' => str_replace($docs_md_base, '', $name),
        ];
    }
}

$sidebar_nav = '<nav class="sidebar-nav"><ul class="ps-0 d-flex flex-column">';
$current_level = 1;
$previous_parent = explode('/', $sidebar_nav_elements[0]['url'])[1];
$nested = false;
foreach($sidebar_nav_elements as $key => $nav) {
    $level = substr_count($nav['url'], '/');
    $path = explode('/', $nav['url']);
    $parent = array_reverse($path)[1];
    $active = $md_fn == $nav['url'] ? 'active' : '';
    while ($level > $current_level) {
        global $nested;
        $nested = true;
        $current_level += 1;
        $show = strpos($md_fn, implode('/',array_slice($path,0,-1))) !== false ? 'show' : '';
        $is_open = $show == 'show' ? 'true': 'false';
        $id = str_replace(" ", "-", strtolower($nav['title']));
        $id = str_replace(":", "", $id);
        $sidebar_nav .= '<button class="display-1 btn d-inline-flex align-items-center rounded " data-bs-toggle="collapse" data-bs-target="#' . $id . '" aria-expanded="'.$is_open.'" aria-current="'.$is_open.'">
                            <i class="fas fa-angle-right me-3"></i><strong>'. ucwords($nav['parent']) .'
                        </strong></button>';
        $sidebar_nav .= '<nav class="collapse ' . $show . '" id="' . $id . '"><ul class="list-unstyled fw-normal pb-1 ps-3 small">';
    }

    if($level == $current_level && $parent != $previous_parent && !$nested) {
        $previous_parent = $parent;
        $sidebar_nav .= '</ul></nav>';
        $show = strpos($md_fn, $parent) !== false ? 'show' : '';
        $is_open = $show == 'show' ? 'true' : 'false';
        $id = str_replace(" ", "-", strtolower($nav['title']));
        $id = str_replace(":", "", $id);
        $sidebar_nav .= '<button class="btn d-inline-flex align-items-center rounded" data-bs-toggle="collapse" data-bs-target="#' . $id . '" aria-expanded="'.$is_open.'" aria-current="'.$is_open. '">
                            <i class="fas fa-angle-right me-3"></i><strong>' . ucwords($nav['parent']) . '
                        </strong></button>';
        $sidebar_nav .= '<nav class="collapse ' . $show . '" id="' . $id . '"><ul class="list-unstyled fw-normal pb-1 ps-3 small">';
    }
    
    
    $sidebar_nav .= '<li class="mb-2"><a href="/' . $nav['url'] . '" class=" '. $active . '">' . $nav['title'] . '</a></li>';
    while ($level < $current_level && $current_level > 1) {
        $nested = false;
        $current_level = $level;
        $sidebar_nav .= '</ul></nav>';
    }
    
}
$sidebar_nav .= '</ul></nav>';

include('../includes/header.php');

?>
<div class="container-xxl main-content">
<?php
    $main_content = '<div class="row flex-wrap-reverse flex-lg-wrap">';

    # right sidebar
    $main_content .= '<div class="col-12 col-lg-3">';
    $main_content .=  $sidebar_nav;
    $main_content .= '</div>';
    $main_content .= '<div class="col-12 col-lg-6"><div class="rendered-markdown">' . $content . '</div>
                </div>';
    # right sidebar
    $main_content .= '<div class="col-12 col-lg-3 ps-2"><div class="side-sub-subnav sticky-top">';
    # ToC
    $main_content .= '<nav class="toc auto-toc pt-2 flex-column border-start">';
    $main_content .= '<strong class="ms-3 d-inline-block w-100 text-secondary border-bottom">On this page</strong>';
    $main_content .= generate_toc($content);
    $main_content .=  '<p class="small text-end mt-3"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
    $main_content .=  '</nav>';

    $main_content .= '</div></div>'; # end of the sidebar col
    $main_content .=  '</div>'; # end of the row

    echo $main_content;
?>
</div>

<?php
include('../includes/footer.php');
