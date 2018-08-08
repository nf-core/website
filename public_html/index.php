<?php include('../includes/header.php'); ?>


    <div class="homepage-header">
      <div class="homepage-header-contents col-md-5">
        <h1><img src="assets/img/logo/nf-core-logo.svg" /></h1>
        <p class="lead font-weight-normal">A community effort to collect a curated set of analysis pipelines built using Nextflow.</p>
        <div class="hompage-cta-flex">
          <a class="hompage-cta" href="/pipelines">View Pipelines</a>
        </div>
      </div>
    </div>

    <div class="triangle triangle-up"></div>
    <div class="homepage-intro">
      <div class="container">
        <div class="row">
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <h3 class="display-5 text-white">For facilities</h3>
            <p class="lead">Highly optimised pipelines with excellent reporting. Validated releases ensure reproducibility.</p>
          </div>
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <h3 class="display-5 text-white">For users</h3>
            <p class="lead">Portable, documented and easy&nbsp;to&nbsp;use workflows.<br>Pipelines that you can trust.</p>
          </div>
          <div class="col-md-4 py-3 py-md-5 px-md-4">
            <h3 class="display-5 text-white">For developers</h3>
            <p class="lead">Companion templates and tools help to validate your code and simplify common tasks.</p>
          </div>
        </div>
      </div>
    </div>
    <div class="triangle triangle-down"></div>

    <div class="container py-5 text-center lead">
      <p>Nextflow is an incredibly powerful and flexible workflow language.</p>
      <p><strong>nf-core</strong> pipelines adhere to strict guidelines - if one works, they all will.</p>
    </div>

    <div class="container homepage-feature-boxes">
      <div class="row">
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Documentation</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/docs.svg" />
              <p class="card-text">Extensive documentation covering installation, usage and description of output files
              ensures that you won't be left in the dark.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">CI Testing</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/travis-ci.svg" />
              <p class="card-text">Every time a change is made to the pipeline code,
              nf-core pipelines use continuous-integration testing to ensure that nothing has broken.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Stable Releases</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/releases.svg" />
              <p class="card-text">nf-core pipelines use GitHub releases to tag stable versions of the code
              and software, making pipeline runs totally reproducable.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title text-center">Docker</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/docker.svg" />
              <p class="card-text">Software dependencies are always available in a bundled docker container,
              which Nextflow can automatically download from dockerhub.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Singularity</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/singularity.svg" />
              <p class="card-text">If you're not able to use Docker, built-in support for Singularity can
              solve your HPC container problems. These are built from the docker containers.</p>
            </div>
          </div>
        </div>
        <div class="col-md-4 mb-5">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">Bioconda</h5>
              <img class="float-right ml-2" height="100px" src="assets/img/bioconda.svg" />
              <p class="card-text">Where possible, pipelines come with a bioconda environment file,
              allowing you to set up a new environment for the pipeline in a single command.</p>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="bg-dark py-5">
      <div class="container">
        <div class="row">
          <div class="col-sm-6">
            <h2 class="text-white">Get started in minutes</h2>
            <p class="text-white-50">Nextflow lets you run nf-core pipelines on virtually any computing environment.</p>
            <p class="text-white-50">nf-core pipelines come with built-in support for
              <a href="https://ewels.github.io/AWS-iGenomes/" target="_blank">AWS iGenomes</a>
              with common species.</p>
            <p class="text-white-50">The nf-core companion tool makes it easy to list all available nf-core pipelines
              and shows which are available locally. Local versions are checked against the latest available release.</p>
          </div>
          <div class="col-sm-6">
            <div class="card text-white bg-dark border-light">
              <div class="card-body">
                <pre class="text-white mb-0">

<span class="text-white-50"># Prerequisites</span>
<span class="text-white-50">* Java 7 or 8 </span>
<span class="text-white-50">* Docker engine 1.10.x (or higher)</span>

<span class="text-white-50"># Install nextflow</span>
curl -s https://get.nextflow.io | bash

<span class="text-white-50"># move nextflow file to your 'PATH'</span>
<span class="text-white-50"># Launch the RNAseq pipeline</span>
nextflow run nf-core/RNAseq \
    -with-docker \
    --genome GRCh37 \
    --reads "data/*_{R1,R2}.fastq.gz"

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
          <a class="btn btn-lg btn-success arrow-hover" href="/pipelines"><span>Available pipelines</span></a>
        </div>
        <div class="col-md-4 px-md-5 mb-5 mb-md-0">
          <h3>Run a pipeline</h3>
          <p>Read the quickstart tutorial to learn how to get set up with the required software and tools, and how to launch a nf-core pipeline.</p>
          <a class="btn btn-lg btn-success arrow-hover" href="/docs"><span>Quickstart Tutorial</span></a>
        </div>
        <div class="col-md-4 px-md-5 mb-5 mb-md-0">
          <h3>Get into the code</h3>
          <p>If you're interested in contributing to nf-core, take a look at the developer documentation to see what's required.</p>
          <a class="btn btn-lg btn-success arrow-hover" href="/docs"><span>Developer docs</span></a>
        </div>
      </div>
    </div>

    <div class="homepage-usedby">
      <div class="container py-5">
        <div class="row">
          <div class="col-12 col-lg-6">
            <h2>Used by groups all over the world</h2>
            <p>The nf-core community is spread all over the globe and includes a large
              number of contributing users. <a href="#">See all &raquo;</a></p>
            <img src="assets/img/institute_images/NGI.svg" data-placement="bottom" data-toggle="tooltip" title="SciLifeLab NGI">
            <img src="assets/img/institute_images/QBiC.svg" data-placement="bottom" data-toggle="tooltip" title="Quantitative Biology Center, Tübingen, Germany">
            <img src="assets/img/institute_images/SciLifeLab.svg" data-placement="bottom" data-toggle="tooltip" title="SciLifeLab, Sweden">
            <img src="assets/img/institute_images/GIS.svg" data-placement="bottom" data-toggle="tooltip" title="GIS, Singapore">
            <img src="assets/img/institute_images/IARC.svg" data-placement="bottom" data-toggle="tooltip" title="IARC, World Health Organisation, France">
          </div>
        </div>
      </div>
    </div>

<?php include('../includes/footer.php'); ?>
