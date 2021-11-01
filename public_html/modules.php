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
            <div class="card 100 pipeline">
                <div class="card-body clearfix">
                    <h3 class="card-title mb-0">
                        <a href="https://github.com/nf-core/modules/tree/master/<?php echo $wf['github_path']; ?>" class="pipeline-name">
                            <span class="d-none d-lg-inline"></span><?php echo $wf['name']; ?>
                        </a> 
                        <button class="btn btn-sm btn-outline-info float-end" data-bs-toggle="collapse" href=".<?php echo $wf['name']; ?>-description" aria-expanded="false">
                            <i class="fas fa-question-circle"></i> descriptions
                        </button>
                    </h3>
                    <p class="card-text mb-1">
                        <?php echo $wf['description']; ?>
                    </p>
                        <?php 
                        $tool_text = "";
                        foreach ($wf['tools'] as $tool) {
                            foreach ($tool as $name => $tool_value) {
                                if($wf['name']!=$name){
                                    $documentation =  $tool_value["documentation"] ? $tool_value["documentation"] : $tool_value["homepage"];
                                    $tool_text .= "<a class=\"ms-2\" data-bs-toggle=\"tooltip\" title=\"documentation\" href=". $documentation ."><i class=\"fas fa-book\"></i></a>";
                                    $tool_text = "<h4 class=\"mb-0\">" . $name . $tool_text. "</h4>";
                                    $description = $tool_value["description"];
                                    $tool_text .= "<span class=\"text-small\" >" .  $description . "</span>";
                                }
                                
                            }
                        }
                        echo $tool_text;
                        ?>
                    <div class="row pt-3">
                        <div class="col-6 border-end border-1">
                            <h4>Input</h4>
                            
                            <?php
                                $input_text = "";
                                foreach ($wf['input'] as $input) {
                                    $input_text .= "<ul class=\"list-switch\">";
                                    foreach ($input as $name => $input_value) {
                                        $input_text .= "<li>";
                                        $input_text .= "<strong>" . $name . " </strong>";
                                        $input_text .= "<span class=\"text-secondary\"> (" . $input_value["type"]. ") </span>";
                                        $description = $input_value["description"];
                                        $description = str_replace("[","<code class=\"px-0\">[",$description);
                                        $description = str_replace("]","]</code>",$description);
                                        $input_text .= "<p class=\"text-small collapse mb-1 " . $wf['name']."-description\" >" .  $description . "</p>";
                                        $input_text .= "</li>";
                                    }
                                    $input_text .= "</ul>";
                                }
                                echo $input_text;
                            ?>
                        </div>
                        <div class="col-6 ">
                            <h4>Output</h4>
                            <?php
                                $output_text = "";
                                foreach ($wf['output'] as $output) {
                                    $output_text .= "<ul class=\"list-switch\">";
                                    foreach ($output as $name => $output_value) {
                                        $output_text .= "<li>";
                                        $output_text .= "<strong>" . $name . " </strong>";
                                        $output_text .= "<span class=\"text-secondary\"> (" . $output_value["type"]. ") </span>";
                                        $description = $output_value["description"];
                                        $description = str_replace("[","<code class=\"px-0\">[",$description);
                                        $description = str_replace("]","]</code>",$description);
                                        $output_text .= "<p class=\"text-small collapse mb-1 " . $wf['name']."-description\" >" .  $description . "</p>";
                                        $output_text .= "</li>";
                                    }
                                    $output_text .= "</ul>";
                                }
                                echo $output_text;
                            ?>
                        </div>
                    </div>
                    </div>
                    <div class="card-footer text-muted">
                        <?php if (count($wf['keywords']) > 0) : ?>
                            <p class="topics mb-0">
                                <?php foreach ($wf['keywords'] as $keyword) : ?>
                                    <a href="/modules?q=<?php echo $keyword; ?>" class="badge pipeline-topic"><?php echo $keyword; ?></a>
                                <?php endforeach; ?>
                                <?php foreach ($wf['authors'] as $author) : ?>
                                    <img class="float-end contrib-avatars" src="https://github.com/<?php echo str_replace("@", "", $author); ?>.png">
                                <?php endforeach; ?>
                            </p>
                        <?php endif; ?>
                    </div>
            </div>
        </div>
    <?php endforeach; ?>
</div>

<?php include('../includes/footer.php'); ?>