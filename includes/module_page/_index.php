<?php

require_once('../includes/functions.php');

########
## Configure page header
########
$title = 'modules/<br class="d-sm-none">' . $module['name'];
$subtitle = $module['description'];
$content = '';
$schema_content = '';
$import_chartjs = true;
$no_auto_toc = true;
$gh_url = "https://github.com/nf-core/modules/tree/master/" . str_replace("/meta.yml", "", $module['github_path']);

# Header - keywords
$header_html = '<p class="mb-0">';

foreach ($module['keywords'] as $keyword) {
    $header_html .= '<a href="/modules?q=' . $keyword . '" class="badge module-topic">' . $keyword . '</a> ';
}
$header_html .= '</p>';

// Highlight any search terms if we have them
if (isset($_GET['q']) && strlen($_GET['q'])) {
    $title = preg_replace("/(" . $_GET['q'] . ")/i", "<mark>$1</mark>", $title);
    $subtitle = preg_replace("/(" . $_GET['q'] . ")/i", "<mark>$1</mark>", $subtitle);
    $header_html = preg_replace("/(" . $_GET['q'] . ")/i", "<mark>$1</mark>", $header_html);
}
# Set defaults (Readme tab)
$pagetab = ''; # empty string is home / readme


# Main page nav and header
$no_print_content = true;
$mainpage_container = false;
include('../includes/header.php');

########
# Extra HTML for the header - tags and GitHub URL
########
?>

<div class="container-fluid mainpage-subheader-heading chevron-down">
    <div class="container text-center">
        <p><a href="<?php echo $gh_url; ?>" class="subheader-link">
                <i class="fab fa-github"></i><?php echo $gh_url; ?>

            </a></p>
    </div>
</div>

<div class="container-xxl main-content">

    <ul class="nav nav-fill nfcore-subnav justify-content-around">
        <li class="nav-item">
            <a class="nav-link<?php if ($pagetab == '') {
                                    echo ' active';
                                } ?>" href="<?php echo $url_base; ?>"><i class="fas fa-sign-in me-1"></i> Introduction</a>
        </li>
    </ul>

    <?php
    ########
    # Make a row with a column for content
    ########
    echo '<div class="row flex-wrap-reverse flex-lg-wrap ms-lg-5"><div class="col-12 col-lg-9">';

    ########
    # Print content
    ########

    ?>
    <div class="module module-page-content mb-2">
        <div class="">
            <div class="module-params ">
                <div class="module-input">
                    <h2 class="text-success">Input</h2>

                    <?php
                    $input_text = '<div class="d-flex flex-column mb-3">';
                    foreach ($module['input'] as $input) {
                        foreach ($input as $name => $input_value) {
                            $description = $input_value['description'];
                            $description = str_replace('[', '<code class="px-0">[', $description);
                            $description = str_replace(']', ']</code>', $description);
                            $input_text .= '<div >';
                            $input_text .= '<span data-bs-toggle="tooltip" title="' . $input_value['description'] . '">' . $name . ' </span>';
                            $input_text .= '<span class="text-muted"> (' . $input_value['type'] . ')</span>';
                            $input_text .= '<p class="text-small collapse show mb-1  ms-3 description ' . $module['name'] . '-description" >' .  $description . '</p>';
                            if (key(end($module['input'])) != $name) { //don't add a comma after the last element
                                $input_text .= '<span class="hide-collapse">,&nbsp;</span>';
                            }
                            $input_text .= '</div>';
                        }
                    }
                    $input_text .= '</div>';
                    echo $input_text;
                    ?>
                </div>
                <div class="module-output">
                    <h2 class="text-success">Output</h2>
                    <?php
                    $output_text = '<div class="d-flex flex-column mb-3">';
                    foreach ($module['output'] as $output) {
                        foreach ($output as $name => $output_value) {
                            $description = $output_value['description'];
                            $description = str_replace('[', '<code class="px-0">[', $description);
                            $description = str_replace(']', ']</code>', $description);
                            $output_text .= '<div >';
                            $output_text .= '<span data-bs-toggle="tooltip" title="' . $output_value['description'] . '">' . $name . ' </span>';
                            $output_text .= '<span class="text-muted"> (' . $output_value['type'] . ')</span>';
                            $output_text .= '<p class="text-small collapse show mb-1 ms-3 description ' . $module['name'] . '-description" >' .  $description . '</p>';
                            if (key(end($module['output'])) != $name) { //don't add a comma after the last element
                                $output_text .= '<span class="hide-collapse">,&nbsp;</span>';
                            }
                            $output_text .= '</div>';
                        }
                    }
                    $output_text .= '</div>';
                    echo $output_text;
                    ?>
                </div>
            </div>
            <div class="module-tools">
                <?php
                $tool_text = '<div class="d-flex flex-wrap">';
                $tool_text .= '<h2 class="text-success">Tools</h2>';
                foreach ($module['tools'] as $tool) {
                    // catch incorrectly formatted yamls
                    if (isset($tool['documentation'])) {
                        $tmp_tool = [];
                        $tmp_tool[key($tool)] = array_slice($tool, 1);
                        $tool = $tmp_tool;
                    }
                    foreach ($tool as $name => $tool_value) {
                        $tool_text .= '<div>';
                        $documentation =  $tool_value['documentation'] ? $tool_value['documentation'] : $tool_value['homepage'];
                        $documentation = '<a class="ms-2" data-bs-toggle="tooltip" title="documentation" href=' . $documentation . '><i class="fas fa-book"></i></a>';
                        $tool_text .= '<h4>' . $name . $documentation . ' </h4>';
                        $description = $tool_value['description'];
                        $tool_text .= '<span class="text-small collapse show description ' . $module['name'] . '-description" >' .  $description . '</span>';
                        $tool_text .= '</div>';

                        $tool_text .= '</ul>';
                    }
                }
                echo $tool_text;
                ?>
        </div>
    </div>
</div>
</div>

<?php
echo '</div>'; # end of the content div
echo '<div class="col-12 col-lg-3 ps-2"><div class="side-sub-subnav sticky-top">';

# module homepage & releases - key stats
if (in_array($pagetab, [''])) {
    // require_once('sidebar.php');
?>
    <div class="module-sidebar">
        <div class="row border-bottom pb-2">
            <div class="col-12">
                <h6><i class="fas fa-terminal fa-xs"></i> command</h6>
                <div class="input-group input-group-sm module-install-cmd">
                    <input type="text" class="form-control input-sm code" id="module-install-cmd-text" data-autoselect="" value="nf-core modules install <?php echo $module['name']; ?>" aria-h3="Copy install command" readonly="">
                    <button class="btn btn-outline-secondary copy-txt" data-bs-target="module-install-cmd-text" data-bs-toggle="tooltip" data-bs-placement="left" title="Copy to clipboard" type="button"><i class="fas fa-clipboard px-1"></i></button>
                </div>
            </div>
        </div>

        <div class="row border-bottom">
            <div class="col-12 contrib-avatars">
                <h6>collaborators</h6>
                <?php foreach ($module['authors'] as $author) : ?>
                    <img src="https://github.com/<?php echo trim(str_replace("@", "", $author)); ?>.png">
                <?php endforeach; ?>
            </div>
        </div>
        <div>
            <h6>get in touch</h6>
            <p><a class="btn btn-sm btn-outline-info" href="https://nfcore.slack.com/channels/modules"><i class="fab fa-slack me-1"></i> Ask a question on Slack</a></p>
            <p><a class="btn btn-sm btn-outline-secondary" href="https://github.com/nf-core/modules/issues"><i class="fab fa-github me-1"></i> Open an issue on GitHub</a></p>
        </div>
    </div>
<?php
}
# Documentation - ToC

echo '</div></div>'; # end of the sidebar col
echo '</div>'; # end of the row

include('../includes/footer.php');
