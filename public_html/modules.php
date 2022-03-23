<?php
require_once('../includes/parse_md.php');

$config = parse_ini_file("../config.ini");
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// get all modules
$sql = "SELECT * FROM nfcore_modules ORDER BY LOWER(name)";
$modules = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        while ($row = mysqli_fetch_array($result)) {
            $row['keywords'] = explode(';', $row['keywords']);
            $row['authors'] = explode(';', $row['authors']);
            $row['tools'] = json_decode($row['tools'], true);
            $row['input'] = json_decode($row['input'], true);
            $row['output'] = json_decode($row['output'], true);
            $modules[] = $row;
        }
        // Free result set
        mysqli_free_result($result);
    } else {
        echo "Oops! Something went wrong. Please try again later.";
    }
}
// get all keywords
$sql = "SELECT keywords FROM nfcore_modules ";
$keywords = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        while ($row = mysqli_fetch_array($result)) {
            $keywords[] = explode(';', $row['keywords']);
        }
        // Free result set
        mysqli_free_result($result);
    } else {
        echo "Oops! Something went wrong. Please try again later.";
    }
}
$keywords_tmp = [];
// flatten keywords array
array_walk_recursive($keywords, function ($v) use (&$keywords_tmp) {
    $keywords_tmp[] = $v;
});
$keywords = $keywords_tmp;

// Close connection
mysqli_close($conn);


// Pagination
// $num_elements = 30;
$num_elements = 3000; // TODO: remove after we have decided on a pagination strategy
$total_pages = ceil(count($modules) / $num_elements);
if (!isset($_GET['page'])) {
    $current_page = 1;
} else {
    $current_page = $_GET['page'];
}
$current_element = ($current_page - 1) * $num_elements;
$current_modules = array_slice($modules, $current_element, $num_elements);

function add_update_url_param($param_key, $param_value)
{
    $params = array_merge($_GET, array($param_key => $param_value));
    return http_build_query($params);
}

$msg = '';
if (isset($_GET['update'])) {
    $output = shell_exec("php ../update_module_details.php 2>&1 | tee -a /home/nfcore/update.log");
    $msg = '<div class="alert alert-success">Manual modules update sync triggered: <pre>' . $output . '</pre></div>';
}


$title = 'Modules';
$subtitle = 'Browse the <strong>' . count($modules) . '</strong> modules that are currently available as part of nf-core.';
include('../includes/header.php');
echo $msg;
?>

<h1>Available Modules</h1>
<p class="mb-3">Modules are the building stones of all DSL2 nf-core blocks. You can find more info <a href="/developers/modules"></a>, if you would like to write your own module.
</p>
<div class=" btn-toolbar mb-4 modules-toolbar" role="toolbar">
    <div class="module-filters input-group input-group-sm w-25">
        <input type="search" class="form-control" placeholder="Search modules" value="<?php echo isset($_GET['q']) ? $_GET['q'] : ''; ?>">
    </div>
</div>
<div class="row flex-wrap-reverse flex-lg-wrap me-lg-5">
    <div class="col-12 col-lg-3 pe-2">
        <div class="facet-bar">
            <?php $keywords_value = array_count_values($keywords);
            arsort($keywords_value);
            ?>
            <ul class="list-unstyled">
                <?php foreach ($keywords_value as $idx => $keyword) : ?>
                    <li class="facet-item" id="<?php echo 'keyword-' . preg_replace('/\s+/', '__', trim($idx)); ?>">
                        <span class=" facet-name"><?php echo trim($idx); ?></span>
                        <span class="facet-value badge rounded-pill bg-secondary float-end">
                            <?php echo $keyword; ?>
                        </span>
                    </li>
                <?php endforeach; ?>
            </ul>
        </div>
        <div class="pipeline_list">
            <ul class="list-unstyled">
                <?php foreach ($pipelines as $pipeline) : ?>
                    <li class="facet-item">
                        <span class="facet-name"><?php echo trim($idx); ?></span>
                        <span class="facet-value badge rounded-pill bg-secondary float-end">
                            <?php echo $keyword; ?>
                        </span>
                    </li>
                <?php endforeach; ?>
            </ul>
        </div>
    </div>
    <div class="col-12 col-lg-9">
        <p class="no-modules text-muted mt-5" style="display: none;">No modules found..</p>
        <div class="modules-container modules-container-list">
            <?php foreach ($current_modules as $idx => $module) : ?>
                <div class="card pipeline module mb-2">
                    <div class="card-body">
                        <div class="module-name">
                            <h4 class="card-title mb-0" id="<?php echo $module['name']; ?>">
                                <a href="modules/<?php echo $module['name']; ?>" class="pipeline-name">
                                    <?php echo $module['name']; ?>
                                </a>

                            </h4>
                        </div>
                        <?php if (count($module['keywords']) > 0) : ?>
                            <div class="topics mb-0">
                                <?php foreach ($module['keywords'] as $keyword) : ?>
                                    <a class="badge  pipeline-topic" data-keyword="<?php echo 'keyword-' . preg_replace('/\s+/', '__', trim($keyword)); ?>"><?php echo $keyword; ?></a>
                                <?php endforeach; ?>
                            </div>
                        <?php endif; ?>
                        <p class="card-text mb-0 mt-2 description <?php echo $module['name']; ?>-description">
                            <?php echo parse_md($module['description'])['content']; ?>
                        </p>
                        <div class="module-params d-flex">
                            <div class="module-input flex-fill">
                                <label class="text-success">Input</label>

                                <?php
                                $input_text = '<p class="d-flex flex-wrap">';
                                foreach ($module['input'] as $input) {
                                    foreach ($input as $name => $input_value) {
                                        $description = $input_value['description'];
                                        $description = str_replace('[', '<code class="px-0">[', $description);
                                        $description = str_replace(']', ']</code>', $description);
                                        $input_text .= '<code class="me-1 mt-1 py-0" data-bs-toggle="tooltip" title="' . $input_value['description'] . '">' . $name . '</code>';
                                    }
                                }
                                $input_text .= '</p>';
                                echo $input_text;
                                ?>
                            </div>
                            <div class="module-output flex-fill">
                                <label class="text-success">Output</label>
                                <?php
                                $output_text = '<p class="d-flex flex-wrap">';
                                foreach ($module['output'] as $output) {
                                    foreach ($output as $name => $output_value) {
                                        $description = $output_value['description'];
                                        $output_text .= '<code class="me-1 mt-1 py-0" data-bs-toggle="tooltip" title="' . $output_value['description'] . '">' . $name . ' </code>';
                                    }
                                }
                                $output_text .= '</p>';
                                echo $output_text;
                                ?>
                            </div>
                            <div class="module-tools  flex-fill">
                                <?php
                                $tool_text = '<div class="d-flex flex-wrap">';
                                foreach ($module['tools'] as $tool) {
                                    // catch incorrectly formatted yamls
                                    if (isset($tool['documentation'])) {
                                        $tmp_tool = [];
                                        $tmp_tool[key($tool)] = array_slice($tool, 1);
                                        $tool = $tmp_tool;
                                    }
                                    foreach ($tool as $name => $tool_value) {
                                        // don't print tools if it has the same name (and therefore usually same description) as the module
                                        if ($module['name'] != $name) {

                                            $tool_text .= '<div>';
                                            $tool_text .= '<span>' . $name . ': </span>';
                                            $description = $tool_value['description'];
                                            $description = parse_md($description)['content'];
                                            $tool_text .= '<span class="text-small description ' . $module['name'] . '-description" >' .  $description . '</span>';
                                            $tool_text .= '</div>';
                                        }
                                        $tool_text .= '</ul>';
                                    }
                                }
                                if (substr_count($tool_text, '<div>') > 0) {
                                    $tool_text = '<label class="text-success">Tools</label>' . $tool_text;
                                }
                                echo $tool_text;
                                ?>
                            </div>
                        </div>
                    </div>
                </div>
        </div>
    <?php endforeach; ?>
    </div>
    <nav aria-label="Module page navigation ">
        <ul class="pagination">
            <?php
            $disable = $current_page - 1 == 0 ? 'disabled' : '';
            $params = add_update_url_param('page', $current_page - 1);
            echo '<li class="page-item ' . $disable . '"><a class="page-link" href= "/modules?' . $params . '">Previous</a></li>';
            for ($page_number = 1; $page_number <= $total_pages; $page_number++) {
                $params = add_update_url_param('page', $page_number);
                $active = $current_page == $page_number ? 'active' : '';
                echo '<li class="page-item ' . $active . '"><a class="page-link" href= "/modules?' . $params . '">' . $page_number . ' </a></li>';
            }
            $params = add_update_url_param('page', $current_page + 1);
            $disable = $current_page - $total_pages == 0 ? 'disabled' : '';
            echo '<li class="page-item ' . $disable . '"><a class="page-link" href= "/modules?' . $params . '">Next</a></li>';
            ?>
        </ul>
    </nav>
</div>
<p class="mt-5 small text-muted">
    <a href="/modules?update">Click here</a> to trigger an update.
</p>
<?php include('../includes/footer.php');
