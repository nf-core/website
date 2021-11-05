<?php
$config = parse_ini_file("../config.ini");
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// Attempt select query execution
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

// Close connection
mysqli_close($conn);

$title = 'Modules';
$subtitle = 'Browse the <strong>' . count($modules) . '</strong> modules that are currently available as part of nf-core.';
include('../includes/header.php');

?>

<h1>Available Modules</h1>
<div class="btn-toolbar mb-4 pipelines-toolbar" role="toolbar">
    <div class="pipeline-filters input-group input-group-sm ms-2 mt-2">
        <input type="search" class="form-control" placeholder="Search keywords" value="<?php echo isset($_GET['q']) ? $_GET['q'] : ''; ?>">
    </div>

</div>

<p class="no-modules text-muted mt-5" style="display: none;">No modules found..</p>

<div class="pipelines-container modules-container-list">
    <?php foreach ($modules as $idx => $wf) : ?>
        <div class="card module mb-2">
            <div class="card-header module-name">
                <h4 class="card-title mb-0" id="<?php echo $wf['name']; ?>">
                    <a href="https://github.com/nf-core/modules/tree/master/<?php echo $wf['github_path']; ?>" class="pipeline-name">
                        <?php echo $wf['name']; ?>
                        <a href="#<?php echo $wf['name']; ?>" class="header-link"><span class=" fas fa-link"></span></a>
                    </a>
                    <button class="btn btn-sm btn-outline-info float-end" data-bs-toggle="collapse" href=".<?php echo $wf['name']; ?>-description" aria-expanded="false">
                        <i class="fas fa-question-circle"></i> descriptions
                    </button>
                </h4>
            </div>
            <div class="card-body">
                <span class="card-text ms-2 collapse <?php echo $wf['name']; ?>-description">
                    <?php echo $wf['description']; ?>
                </span>
                <div class="module-params d-flex justify-content-between">
                    <div class="module-input">
                        <label class="text-success">Input</label>

                        <?php
                        $input_text = '<div class="d-flex  flex-wrap">';
                        foreach ($wf['input'] as $input) {
                            foreach ($input as $name => $input_value) {
                                $input_text .= '<div >';
                                $input_text .= '<span>' . $name . ' </span>';
                                $input_text .= '<span class="text-muted"> (' . $input_value['type'] . ')</span>';
                                $description = $input_value['description'];
                                $description = str_replace('[', '<code class="px-0">[', $description);
                                $description = str_replace(']', ']</code>', $description);
                                $input_text .= '<p class="text-small collapse mb-1 ' . $wf['name'] . '-description" >' .  $description . '</p>';
                                if (key(end($wf['input'])) != $name) { //don't add a comma after the last element
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
                        <label class="text-success">Output</label>
                        <?php
                        $output_text = '<div class="d-flex  flex-wrap">';
                        foreach ($wf['output'] as $output) {
                            foreach ($output as $name => $output_value) {
                                $output_text .= '<div >';
                                $output_text .= '<span>' . $name . ' </span>';
                                $output_text .= '<span class="text-muted"> (' . $output_value['type'] . ')</span>';
                                $description = $output_value['description'];
                                $description = str_replace('[', '<code class="px-0">[', $description);
                                $description = str_replace(']', ']</code>', $description);
                                $output_text .= '<p class="text-small collapse mb-1 ' . $wf['name'] . '-description" >' .  $description . '</p>';
                                if (key(end($wf['output'])) != $name) { //don't add a comma after the last element
                                    $output_text .= '<span class="hide-collapse">,&nbsp;</span>';
                                }
                                $output_text .= '</div>';
                            }
                        }
                        $output_text .= '</div>';
                        echo $output_text;
                        ?>
                    </div>
                    <div class="module-tools">
                        <?php
                        $tool_text = '<div class="d-flex flex-wrap">';
                        foreach ($wf['tools'] as $tool) {
                            foreach ($tool as $name => $tool_value) {
                                // don't print tools if it has the same name (and therefore usually same description) as the module
                                if ($wf['name'] != $name) {
                                    // catch incorrectly formatted yamls
                                    $tool_text .= '<div>';
                                    $documentation =  $tool_value['documentation'] ? $tool_value['documentation'] : $tool_value['homepage'];
                                    $documentation = '<a class="ms-2" data-bs-toggle="tooltip" title="documentation" href=' . $documentation . '><i class="fas fa-book"></i></a>';
                                    $tool_text .= '<span>' . $name . $documentation . ' </span>';
                                    $description = $tool_value['description'];
                                    $tool_text .= '<span class="text-small collapse ' . $wf['name'] . '-description" >' .  $description . '</span>';
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
        <div class="card-footer text-muted d-flex align-items-center justify-content-between">
            <?php if (count($wf['keywords']) > 0) : ?>
                <div class="topics mb-0">
                    <?php foreach ($wf['keywords'] as $keyword) : ?>
                        <a href="/modules?q=<?php echo $keyword; ?>" class="badge pipeline-topic"><?php echo $keyword; ?></a>
                    <?php endforeach; ?>
                </div>
            <?php endif; ?>
            <div>
                <?php foreach ($wf['authors'] as $author) : ?>
                    <img class="float-end contrib-avatars" src="https://github.com/<?php echo trim(str_replace("@", "", $author)); ?>.png">
                <?php endforeach; ?>
            </div>
        </div>
</div>
<?php endforeach; ?>
</div>

<?php include('../includes/footer.php'); ?>