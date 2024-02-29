<?php
require_once '../includes/functions.php';
require_once '../includes/parse_md.php';

$config = parse_ini_file('../config.ini');
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// get pipelines which use this module
$sql =
    "SELECT DISTINCT nfcore_pipelines.name,nfcore_pipelines.html_url FROM pipelines_modules
        INNER JOIN nfcore_modules ON pipelines_modules.module_id = nfcore_modules.id
        INNER JOIN nfcore_pipelines ON pipelines_modules.pipeline_id = nfcore_pipelines.id
        WHERE nfcore_modules.name = '" .
    $module['name'] .
    "' ORDER BY LOWER(nfcore_pipelines.name)";
$pipelines = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipelines = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}
mysqli_close($conn);
// function to create bootstrap row with the first column for the name and type and the second for the description
function create_row($name, $type, $description, $pattern) {
    $id = strtolower(preg_replace('/[^\w\-\.]+/', '', str_replace(' ', '-', $name)));
    $id = str_replace('.', '-', $id); // periods break the js code, because they are not valid in selector ids

    $row = '<div class="row border-bottom align-items-center">';
    $row .= '<div class="col-12 col-md-3 small-h">';
    $row .= '<h4 id="' . $id . '" class="module-row-name-id">';
    $row .= '<code>' . $name . '</code>';
    $row .= '<span class="text-muted"> (' . $type . ')</span>';
    $row .=
        '<a href=#' .
        $id .
        ' class="header-link scroll_to_link me-2"><span class="fas fa-link" aria-hidden="true"></span></a>';
    $row .= '</h4>';
    $row .= '</div>';
    $row .= '<div class=" col-12 col-md' . ($pattern != '' ? '-5' : '-7') . '">';
    $row .= '<span class="small">' . parse_md($description)['content'] . '</span>';
    $row .= '</div>';
    $row .= '<div class="col-12 col-md' . ($pattern != '' ? '-4' : '-1') . ' ms-auto">';
    if ($pattern != '') {
        $row .= '<code class="float-end">' . $pattern . '</code>';
    }
    $row .= '</div>';

    $row .= '</div>';
    return $row;
}
$header = '<div class="row border-bottom border-3">';
$header .= '<div class="col-12 col-md-3">';
$header .= '<span class="text-muted">Name</span>';
$header .= '</div>';
$header .= '<div class=" col-12 col-md-6">';
$header .= '<span class="text-muted">Description</span>';
$header .= '</div>';
$header .= '<div class="col-12 col-md">';
$header .= '<span class="text-muted float-end">Pattern</span>';
$header .= '</div>';
$header .= '</div>';

########
## Configure page header
########
$title = 'modules/<br class="d-sm-none">' . $module['name'];
$subtitle = parse_md($module['description'])['content'];
$content = '';
$schema_content = '';
$import_chartjs = true;
$no_auto_toc = true;
$gh_url = 'https://github.com/nf-core/modules/tree/master/' . str_replace('/meta.yml', '', $module['github_path']);

# Header - keywords
$header_html = '<p class="mb-0">';

foreach ($module['keywords'] as $keyword) {
    $header_html .= '<a href="/modules?q=' . $keyword . '" class="badge module-topic">' . $keyword . '</a> ';
}
$header_html .= '</p>';

// Highlight any search terms if we have them
if (isset($_GET['q']) && strlen($_GET['q'])) {
    $quoted_search_term = preg_quote($_GET['q'], '/');
    $title = preg_replace('/(' . $quoted_search_term . ')/i', "<mark>$1</mark>", $title);
    $subtitle = preg_replace('/(' . $quoted_search_term . ')/i', "<mark>$1</mark>", $subtitle);
    $header_html = preg_replace('/(' . $quoted_search_term . ')/i', "<mark>$1</mark>", $header_html);
}
# Set defaults (Readme tab)
$pagetab = ''; # empty string is home / readme
# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
include '../includes/header.php';

########
# Extra HTML for the header - tags and GitHub URL
########
?>

<div class="container-fluid mainpage-subheader-heading chevron-down pt-5">
    <div class="container d-flex flex-column align-items-center">
        <p>
        <div class="input-group input-group module-install-cmd w-50">
            <span class="input-group-text"><i class="fas fa-terminal"></i></span>
            <input type="text" class="form-control input code" id="module-install-cmd-text" data-autoselect=""
            value="nf-core modules install <?php echo str_replace('_', '/', $module['name']); ?>"
            aria-h3="Copy install command" readonly="">
            <button class="btn btn-outline-secondary copy-txt" data-bs-target="module-install-cmd-text" data-bs-toggle="tooltip" data-bs-placement="left" title="Copy to clipboard" type="button"><i class="fas fa-clipboard px-1"></i></button>
        </div>
        </p>
        <p><a href="<?php echo $gh_url; ?>" class="subheader-link">
                <i class="fab fa-github"></i> <?php echo $gh_url; ?>
            </a></p>
    </div>
</div>

<div class="container-xxl main-content">

    <!-- <ul class="nav nav-fill nfcore-subnav justify-content-around">
        <li class="nav-item">
            <a class="nav-link<?php if ($pagetab == '') {
                echo ' active';
            } ?>" href="<?php echo $url_base; ?>"><i class="fas fa-sign-in me-1"></i> Introduction</a>
        </li>
    </ul> -->

    <?php ########

# Make a row with a column for content
    ########
    echo '<div class="row flex-wrap-reverse flex-lg-wrap ms-lg-5"><div class="col-12 col-lg-9">';
########
# Print content
########
?>
    <div class="module module-page-content mb-2">
        <h2 id="description" class="ms-n3"><i class=" far fa-book fa-fw"></i> Description
            <a href="#description" class="header-link scroll_to_link"><span class="fas fa-link" aria-hidden="true"></span></a>
        </h2>
        <p class="ps-2"><?php echo parse_md($module['description'])['content']; ?></p>
        <div class="module-params ">
            <div class="module mt-5-input">
                <h2 id="input" class="text-success ms-n3"><i class="fad fa-sign-in fa-fw"></i> Input
                    <a href="#input" class="header-link scroll_to_link"><span class="fas fa-link" aria-hidden="true"></span></a>
                </h2>


                <?php
                $input_text = '<div class="">';
                $input_text .= $header;
                foreach ($module['input'] as $input) {
                    foreach ($input as $name => $input_value) {
                        $description = $input_value['description'];
                        $input_text .= create_row($name, $input_value['type'], $description, $input_value['pattern']);
                    }
                }
                $input_text .= '</div>';
                echo $input_text;
                ?>
            </div>
            <div class="module-output mt-5">
                <h2 id="output" class="text-success ms-n3"><i class="fad fa-sign-out fa-fw"></i> Output
                    <a href="#output" class="header-link scroll_to_link"><span class="fas fa-link" aria-hidden="true"></span></a>
                </h2>

                <?php
                $output_text = '<div class="">';
                $output_text .= $header;
                foreach ($module['output'] as $output) {
                    foreach ($output as $name => $output_value) {
                        $description = $output_value['description'];
                        $output_text .= create_row(
                            $name,
                            $output_value['type'],
                            $description,
                            $output_value['pattern'],
                        );
                    }
                }
                $output_text .= '</div>';
                echo $output_text;
                ?>
            </div>
            <div class="module-tools mt-5">
                <?php
                $tool_text = '<div class="">';
                $tool_text .= '<h2 id="tools" class="text-success ms-n3"><i class="far fa-wrench fa-fw"></i> Tools';
                $tool_text .=
                    '<a href="#tools" class="header-link scroll_to_link"><span class="fas fa-link" aria-hidden="true"></span></a>';
                $tool_text .= '</h2>';
                foreach ($module['tools'] as $tool) {
                    // catch incorrectly formatted yamls
                    if (isset($tool['documentation'])) {
                        $tmp_tool = [];
                        $tmp_tool[key($tool)] = array_slice($tool, 1);
                        $tool = $tmp_tool;
                    }
                    foreach ($tool as $name => $tool_value) {
                        $tool_text .= '<div>';
                        $documentation_url = $tool_value['documentation']
                            ? $tool_value['documentation']
                            : $tool_value['homepage'];
                        $documentation =
                            '<a class="btn btn-outline-secondary float-end" data-bs-toggle="tooltip" title="documentation" href=' .
                            $documentation_url .
                            '><i class="far fa-books"></i> Documentation</a>';
                        $tool_text .= '<h4>' . $name . $documentation . ' </h4>';
                        $description = $tool_value['description'];
                        $description = parse_md($description)['content'];
                        $tool_text .=
                            '<span class="small collapse show description ' .
                            $module['name'] .
                            '-description" >' .
                            $description .
                            '</span>';
                        $tool_text .= '<div class="d-flex mt-3">';
                        if (
                            ($tool_value['tool_dev_url'] != $documentation_url) &
                            ($tool_value['tool_dev_url'] != '') &
                            (trim(strtolower($tool_value['tool_dev_url'])) != 'none')
                        ) {
                            $tool_dev_url = $tool_value['tool_dev_url'];
                            $tool_dev_icon = preg_match('/github.com/', $tool_dev_url)
                                ? 'fab fa-github'
                                : 'far fa-file-code';
                            $tool_dev_icon = '<i class="' . $tool_dev_icon . ' me-1"></i>';
                            $tool_text .=
                                '<span class="badge bg-secondary me-3">' .
                                $tool_dev_icon .
                                '<a class="text-white" href="' .
                                $tool_dev_url .
                                '">' .
                                $tool_dev_url .
                                '</a></span>';
                        }

                        if ($tool_value['doi'] != '') {
                            $tool_text .=
                                '<a class="badge bg-secondary text-white me-3" href="https://doi.org/' .
                                $tool_value['doi'] .
                                '">doi: ' .
                                $tool_value['doi'] .
                                '</a>';
                        }
                        if ($tool_value['licence'] != '') {
                            $tool_text .=
                                '<span class="badge bg-secondary text-white">License: ' .
                                implode(', ', $tool_value['licence']) .
                                '</span>';
                        }
                        $tool_text .= '</div>';
                        $tool_text .= '</div>';
                        $tool_text .= '</ul>';
                    }
                }
                $tool_text .= '</div>';
                echo $tool_text;
                ?>
            </div>
        </div>
    </div>

    <?php
    echo '</div>'; # end of the content div
    echo '<div class="col-12 col-lg-3 ps-3 h-100"><div class="side-sub-subnav">';
    # module homepage & releases - key stats
    if (in_array($pagetab, [''])) { ?>
        <div class="module-sidebar">
            <?php if (count($pipelines) > 0): ?>
                <div class="row border-bottom">
                    <div class="col-12">
                        <h6>nf-core pipelines with this module</h6>
                        <?php foreach ($pipelines as $pipeline) {
                            echo '<a class="badge bg-light text-dark me-1 mb-1" href="/' .
                                $pipeline['name'] .
                                '">' .
                                $pipeline['name'] .
                                '</a>';
                        } ?>
                    </div>
                </div>
            <?php endif; ?>
            <div class="row border-bottom">
                <div class="col-12 contrib-avatars">
                    <h6>Authors</h6>
                    <?php foreach ($module['authors'] as $author): ?>
                        <a href="https://github.com/<?php echo trim(str_replace('@', '', $author)); ?>">
                            <img src="https://github.com/<?php echo trim(str_replace('@', '', $author)); ?>.png">
                        </a>
                    <?php endforeach; ?>
                </div>
            </div>
            <div>
                <h6>get in touch</h6>
                <p><a class="btn btn-sm btn-outline-info" href="https://nfcore.slack.com/channels/modules"><i class="fab fa-slack me-1"></i> Ask a question on Slack</a></p>
                <p><a class="btn btn-sm btn-outline-secondary" href="https://github.com/nf-core/modules/issues"><i class="fab fa-github me-1"></i> Open an issue on GitHub</a></p>
            </div>
        </div>
    <?php }
    echo '</div></div>'; # end of the sidebar col
    echo '</div>';

# end of the row
?>
    <div class="toast cmd_copied" id="module_sidebar_cmd_copied" role="alert" aria-live="assertive" aria-atomic="true">
        <div class="toast-header">
            <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded me-2" alt="">
            <strong class="me-auto">Command copied!</strong>
            <button type="button" class="ms-2 mb-1 btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
        </div>
        <div class="toast-body">
            Paste this command into your terminal to download this module.
        </div>
    </div>
    <?php include '../includes/footer.php';
