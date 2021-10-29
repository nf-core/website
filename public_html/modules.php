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
    <div class="btn-group btn-group-sm mt-2 d-none d-xl-block" role="group">
        <button type="button" class="btn text-body">Display:</button>
    </div>
    <div class="btn-group btn-group-sm mt-2" role="group">
        <button data-dtype="blocks" type="button" class="display-btn btn btn-outline-success active" title="Display as blocks" data-bs-toggle="tooltip"><i class="fas fa-th-large"></i></button>
        <button data-dtype="list" type="button" class="display-btn btn btn-outline-success" title="Display as list" data-bs-toggle="tooltip"><i class="fas fa-bars"></i></button>
    </div>
    <div class="pipeline-filters input-group input-group-sm ms-2 mt-2">
        <input type="search" class="form-control" placeholder="Search keywords" value="<?php echo isset($_GET['q']) ? $_GET['q'] : ''; ?>">
    </div>

</div>

<p class="no-modules text-muted mt-5" style="display: none;">No modules found..</p>

<div class="row row-cols-1 row-cols-md-2 g-4 pipelines-container">

    <?php foreach ($modules as $idx => $wf) : ?>
        <div class="col">
            <div class="card h-100 pipeline">
                <div class="card-body clearfix">
                    <h3 class="card-title mb-0">
                        <a href="https://github.com/nf-core/modules/tree/master/<?php echo $wf['github_path']; ?>" class="pipeline-name">
                            <span class="d-none d-lg-inline"></span><?php echo $wf['name']; ?>
                        </a>
                    </h3>
                    <p class="card-text mb-0 mt-2">
                        <?php echo $wf['description']; ?>
                    </p>
                    <div class="accordion accordion-flush" id="accordion-<?php echo $wf['name']; ?>">
                        <div class="accordion-item">
                            <h2 class="accordion-header text-white" id="heading_tools">
                                <button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#collapse_<?php echo $idx; ?>_tools" aria-expanded="true" aria-controls="collapse_<?php echo $idx; ?>_tools">
                                    tools
                                </button>
                            </h2>
                            <div id="collapse_<?php echo $idx; ?>_tools" class="accordion-collapse collapse" aria-labelledby="heading_tools" data-bs-parent="#accordion-<?php echo $wf['name']; ?>">
                                <div class="accordion-body">
                                    <?php foreach ($wf['tools'] as $tool) {
                                        foreach ($tool as $name => $tool_value) {
                                            if (!is_array($tool_value)) {
                                                $tool_text .= "<dt>" . $name . "</dt>";
                                                $tool_text .= "<dd>" . $tool_value . "</dd>";
                                            }
                                            $tool_text = "<h3>" . $name . "</h3>";
                                            $tool_text .= "<dl>";

                                            foreach ($tool_value as $key => $value) {
                                                if (is_array($value)) {
                                                    $value = implode("; ", $value);
                                                } else if ($value == "None") {
                                                    continue;
                                                }
                                                $tool_text .= "<dt>" . $key . "</dt>";
                                                $tool_text .= "<dd>" . $value . "</dd>";
                                            }
                                            $tool_text .= "</dl>";
                                        }
                                    }
                                    echo $tool_text;
                                    ?>
                                </div>
                            </div>
                        </div>
                        <div class="accordion-item">
                            <h2 class="accordion-header text-white" id="heading_input">
                                <button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#collapse_<?php echo $idx; ?>_input" aria-expanded="true" aria-controls="collapse_<?php echo $idx; ?>_input">
                                    input
                                </button>
                            </h2>
                            <div id="collapse_<?php echo $idx; ?>_input" class="accordion-collapse collapse" aria-labelledby="heading_input" data-bs-parent="#accordion-<?php echo $wf['name']; ?>">
                                <div class="accordion-body">
                                    <?php
                                    $input_text = "";
                                    foreach ($wf['input'] as $input) {
                                        foreach ($input as $name => $input_value) {
                                            $input_text .= "<h4>" . $name . "</h4>";
                                            $input_text .= "<dl>";
                                            foreach ($input_value as $key => $value) {
                                                if (is_array($value)) {
                                                    $value = implode("; ", $value);
                                                } elseif ($key == "pattern") {
                                                    $value = "<code>" . $value . "</code>";
                                                }
                                                $input_text .= "<dt>" . $key . "</dt>";
                                                $input_text .= "<dd>" . $value . "</dd>";
                                            }
                                            $input_text .= "</dl>";
                                        }
                                    }
                                    echo $input_text;
                                    ?>
                                </div>
                            </div>
                        </div>
                        <div class="accordion-item">
                            <h2 class="accordion-header text-white" id="heading_output">
                                <button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#collapse_<?php echo $idx; ?>_output" aria-expanded="true" aria-controls="collapse_<?php echo $idx; ?>_output">
                                    output
                                </button>
                            </h2>
                            <div id="collapse_<?php echo $idx; ?>_output" class="accordion-collapse collapse" aria-labelledby="heading_output" data-bs-parent="#accordion-<?php echo $wf['name']; ?>">
                                <div class="accordion-body">
                                    <?php
                                    $output_text = "";
                                    foreach ($wf['output'] as $output) {
                                        foreach ($output as $name => $output_value) {
                                            $output_text .= "<h4>" . $name . "</h4>";
                                            $output_text .= "<dl>";
                                            foreach ($output_value as $key => $value) {
                                                if (is_array($value)) {
                                                    $value = implode("; ", $value);
                                                } elseif ($key == "pattern") {
                                                    $value = "<code>" . $value . "</code>";
                                                }
                                                $output_text .= "<dt>" . $key . "</dt>";
                                                $output_text .= "<dd>" . $value . "</dd>";
                                            }
                                            $output_text .= "</dl>";
                                        }
                                    }
                                    echo $output_text;
                                    ?>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="card-footer text-muted">
                    <?php if (count($wf['keywords']) > 0) : ?>
                        <p class="topics mb-0">
                            <?php foreach ($wf['keywords'] as $keyword) : ?>
                                <a href="/modules?q=<?php echo $keyword; ?>" class="badge pipeline-topic"><?php echo $keyword; ?></a>
                            <?php endforeach; ?>
                        </p>
                    <?php endif; ?>
                </div>
            </div>
        </div>
    <?php endforeach; ?>
</div>

<?php include('../includes/footer.php'); ?>