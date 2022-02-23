<?php


$title = 'nf-cordle';
$subtitle = 'Guess the pipeline and the modules';
$config = parse_ini_file("../config.ini");
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// get all modules
$sql = "SELECT * FROM nfcore_modules ORDER BY LOWER(name)";
$modules = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        while ($row = mysqli_fetch_array($result)) {
            $modules[] = $row['name'];
        }
        // Free result set
        mysqli_free_result($result);
    } else {
        echo "Oops! Something went wrong. Please try again later.";
    }
}
$pipelines = array(
    "fetchngs" => array(
        "custom_dumpsoftwareversions"
    ),
    "rnaseq" => array(
        "fastqc",
        "cutadapt",
        "fastqc"
    ),
    "test" => array(
        "abacas",
        "abacas"
    )

);
$current_pipeline_name = "test";
$current_pipeline = $pipelines[$current_pipeline_name];
echo '<script type =  "text/javascript" >var current_pipeline_name = "' . $current_pipeline_name . '"; var current_pipeline_modules = "' . implode($current_pipeline, ";") . '";</script>';
include('../includes/header.php');
include('nf-cordle_pipelines.php');
?>
<h1 class="my-5">nf-cordle/modules</h1>
<div class="nf-cordle">
    <div data-guess-grid class="guess-grid">
        <div class="alert-container"></div>
        <?php foreach ($current_pipeline as $pipeline_module) : ?>
            <select class="guess-tile tile-module text-center" placeholder="Select module">
                <option></option>
                <?php foreach ($modules as $module) : ?> <option value="<?php echo $module; ?>"><?php echo $module; ?></option>
                <?php endforeach; ?>
            </select>
        <?php endforeach; ?>
    </div>
    <button class="btn btn-success submit-guess">Submit!</button>

</div>

<?php

include('../includes/footer.php');
