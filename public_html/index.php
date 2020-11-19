<?php
// Get a random subset of contributor institute logos
require_once("../includes/libraries/Spyc.php");
$contributors = spyc_load_file('../nf-core-contributors.yaml');
$contributors_img_list = [];
foreach($contributors['contributors'] as $idx => $c){
  if(isset($c['image_fn']) and $c['image_fn'] and isset($c['full_name']) and $c['full_name']){
    $card_id = preg_replace('/[^a-z]+/', '-', strtolower($c['full_name']));
    $img_path = 'assets/img/contributors-white/'.$c['image_fn'];
    if(file_exists($img_path)){
      $contributors_img_list[] = '<a href="/community#'.$card_id.'"><img src="'.$img_path.'" data-placement="bottom" data-toggle="tooltip" title="'.$c['full_name'].'"></a>';
    }
  }
}
// Shuffle and truncate the list
shuffle($contributors_img_list);

include('../includes/header.php');
?>


    <div class="homepage-header">
      <div class="homepage-header-contents col-md-5">
        <h1>
            <img src="assets/img/logo/nf-core-logo.svg" class="hide-dark">
            <img src="assets/img/logo/nf-core-logo-darkbg.svg" class="hide-light">
        </h1>
        <p class="lead font-weight-normal">A community effort to collect a curated set of analysis pipelines built using Nextflow.</p>
        <div class="hompage-cta-flex">
          <a class="hompage-cta" href="/pipelines">View Pipelines</a>
        </div>
      </div>
    </div>
    <form class="form-inline searchbar_form homepage_search" action="search" method="get">
      <input type="search" class="form-control" placeholder="Search" name="q" required>
      <button type="submit" class="btn btn-outline-success">Search</button>
    </form>

    <div class="triangle triangle-up"></div>
    <div class="homepage-intro">
      <div class="container">
        <div class="row">
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <i class="fad fa-industry-alt fa-3x mb-2"></i>
            <h3 class="display-5 text-white">For facilities</h3>
            <p class="lead">Highly optimised pipelines with excellent reporting. Validated releases ensure reproducibility.</p>
          </div>
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <i class="fad fa-users fa-3x mb-2"></i>
            <h3 class="display-5 text-white">For users</h3>
            <p class="lead">Portable, documented and easy&nbsp;to&nbsp;use workflows.<br>Pipelines that you can trust.</p>
          </div>
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <i class="fad fa-laptop-code fa-3x mb-2"></i>
            <h3 class="display-5 text-white">For developers</h3>
            <p class="lead">Companion templates and tools help to validate your code and simplify common tasks.</p>
          </div>
        </div>
        <p class="text-center text-white m-0 pb-3">
          nf-core is now published in <a class="text-white link-underline" href="https://www.nature.com/articles/s41587-020-0439-x" target="_blank">Nature Biotechnology</a>!
          <a class="text-white link-underline" href="https://rdcu.be/b1GjZ" target="_blank">Read the full text here</a>.
        </p>
      </div>
    </div>
    <div class="triangle triangle-down"></div>

    <div class="container py-5 text-center lead" id="features">
      <p>Nextflow is an incredibly powerful and flexible workflow language.</p>
      <p><strong>nf-core</strong> pipelines adhere to strict guidelines - if one works, they all will.</p>
    </div>

    <div id="features" class="container homepage-feature-boxes">
      <div class="row">
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Documentation</h5>
              <img class="float-right ml-3" height="100px" src="assets/img/docs.svg" />
              <p class="card-text">Extensive documentation covering installation, usage and description of output files
              ensures that you won't be left in the dark.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">CI Testing</h5>
              <img class="float-right ml-3" height="90px" src="assets/img/github-actions.svg" />
              <p class="card-text">Every time a change is made to the pipeline code,
              nf-core pipelines use continuous-integration testing to ensure that nothing has broken.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Stable Releases</h5>
              <img class="float-right ml-3" height="100px" src="assets/img/releases.svg" />
              <p class="card-text">nf-core pipelines use GitHub releases to tag stable versions of the code
              and software, making pipeline runs totally reproducible.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Docker</h5>
              <img class="float-right ml-3" height="100px" src="assets/img/docker.svg" />
              <p class="card-text">Software dependencies are always available in a bundled docker container,
              which Nextflow can automatically download from Docker Hub.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Singularity</h5>
              <img class="float-right ml-3" height="100px" src="assets/img/singularity.svg" />
              <p class="card-text">If you're not able to use Docker, built-in support for Singularity can
              solve your HPC container problems. These are built from the docker containers.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Bioconda</h5>
              <img class="float-right ml-3" height="100px" src="assets/img/bioconda.svg" />
              <p class="card-text">Where possible, pipelines come with a bioconda environment file,
              allowing you to set up a new environment for the pipeline with a single command.</p>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="bg-secondary py-5">
      <div class="container">

        <div class="row videos-row">
          <div class="col-md-6">
            <div class="embed-responsive embed-responsive-16by9 hidden-xs hidden-sm">
              <iframe id="nf-core-video" class="embed-responsive-item" src="https://www.youtube.com/embed/gUM9acK25tQ" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
            </div>
          </div>
          <div class="col-md-6">
            <ul class="list-switch left video-chooser">
              <li><a href="https://youtu.be/gUM9acK25tQ" data-src="https://www.youtube.com/embed/gUM9acK25tQ" class="active"><span class="hidden-lg hidden-md label label-default">Video:</span> Introduction to nf-core <em>(1:01)</em></a></li>
              <!-- <li><a href="https://youtu.be/Gg5neIPuiVo" data-src="https://www.youtube.com/embed/Gg5neIPuiVo"><span class="hidden-lg hidden-md label label-default">Video:</span> Installing MultiQC <em>(4:33)</em></a></li> -->
            </ul>
            <p class="mt-3 ml-md-4 text-white-50 small"><em>More videos coming soon..</em></p>
          </div>
        </div>

      </div>
    </div>

    <div class="bg-dark py-5">
      <div class="container">
        <div class="row">
          <div class="col-sm-6">
            <h2 id="get-started" class="text-white">Get started in minutes</h2>
            <p class="text-white-50">Nextflow lets you run nf-core pipelines on virtually any computing environment.</p>
            <p class="text-white-50">Most nf-core genomics pipelines come with built-in support for
              <a href="https://ewels.github.io/AWS-iGenomes/" target="_blank" style="white-space:nowrap;">AWS-iGenomes</a>,
              with genome references for over 30 common species.</p>
            <p class="text-white-50">The nf-core companion tool makes it easy to list all available nf-core pipelines
              and shows which are available locally. Local versions are checked against the latest available release.</p>
          </div>
          <div class="col-sm-6">
            <div class="card text-white bg-dark border-light">
              <div class="card-body">
                <pre class="text-white mb-0">
<span class="text-white-50"># Install nextflow</span>
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

<span class="text-white-50"># Launch the RNAseq pipeline</span>
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --genome GRCh37 \
    -profile docker

<span class="text-white-50"># Install nf-core tools</span>
pip install nf-core

<span class="text-white-50"># List all nf-core pipelines and show available updates</span>
nf-core list
</pre>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="container my-5 text-center">
      <p class="lead"></p>
      <div class="row">
        <div class="col-md-4 px-md-5 mb-5 mb-md-0">
          <h3>See what's available</h3>
          <p>Check out the available pipelines to see if we have what you need. Each comes with release details, keywords and a description.</p>
          <a class="btn btn-lg btn-success arrow-hover" href="/pipelines"><span><i class="fas fa-toolbox mr-1"></i> Available pipelines</span></a>
        </div>
        <div class="col-md-4 px-md-5 mb-5 mb-md-0">
          <h3>Run a pipeline</h3>
          <p>Read the quickstart tutorial to learn how to get set up with the required software and tools, and how to launch a nf-core pipeline.</p>
          <a class="btn btn-lg btn-success arrow-hover" href="/usage/introduction"><span><i class="fad fa-graduation-cap mr-1"></i> Quickstart Tutorial</span></a>
        </div>
        <div class="col-md-4 px-md-5 mb-5 mb-md-0">
          <h3>Get into the code</h3>
          <p>If you're interested in contributing to nf-core, take a look at the developer documentation to see what's required.</p>
          <a class="btn btn-lg btn-success arrow-hover" href="developers/guidelines"><span><i class="fad fa-code mr-1"></i> Developer docs</span></a>
        </div>
      </div>
    </div>

    <div id="community" class="homepage-usedby">
      <div class="container py-5">
        <h2>
          <a class="btn btn-success float-right d-none d-md-inline" href="/community#organisations">See a complete list &raquo;</a>
          <a href="/community#organisations">Used by groups all over the world</a>
        </h2>
        <p>The nf-core community is spread all over the globe and includes a large
          number of contributing users.</p>
        <p><a class="btn btn-success d-inline d-md-none" href="/community#organisations">See a complete list &raquo;</a></p>
        <div class="homepage_contrib_logos">
          <?php foreach($contributors_img_list as $idx => $img){
            // Hide images after 18 shown
            if($idx > 16) echo str_replace('<a href', '<a style="display:none;" href', $img);
            else echo str_replace('<a href', '<a class="contrib_shown" href', $img);
          } ?>
        </div>
      </div>
    </div>

<?php include('../includes/footer.php'); ?>
