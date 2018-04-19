<?php
$title = 'Documentation';
$subtitle = 'Learn how to develop pipelines that adhere to the nf-core guidelines.';
include('../includes/header.php');
?>


<h1 id="introduction">Introduction</h1>
<p>These pages are still very much under construction and are likely to change a lot in the near future. If you have thoughts please join the discussion!</p>

<h2 id="getting-help">Getting Help</h2>
<p>The quickest place to get help is on the nf-core Gitter channel: <a href="https://gitter.im/nf-core/Lobby">https://gitter.im/nf-core/</a> - a free chat interface that integrates nicely with GitHub.</p>

<h1 id="guidelines-for-new-pipelines">Guidelines for New Pipelines</h1>
<p>For now, the guidelines are pretty simple. The below gives an outline of what is required.
For more detail, please see the <a href="errors"><strong>list of lint test error codes</strong></a>.</p>

<p>If in doubt, we recommend using the <a href="https://github.com/nf-core/cookiecutter">nf-core/cookiecutter</a>
pipeline template to get started. This generates a skeleton pipeline with all of the features
required to pass the below tests.</p>

<h2 id="features">Features</h2>
<p>All pipelines must adhere to the following:</p>

<ul>
  <li>Be built using Nextflow</li>
  <li>An <a href="https://choosealicense.com/licenses/mit/">MIT licence</a></li>
  <li>Software bundled using <a href="https://www.docker.com/">Docker</a>
    <ul>
      <li>This must be at least one <code class="highlighter-rouge">Dockerfile</code> in the repository</li>
      <li>Automatic builds will be set up with tagged versions for GitHub releases</li>
      <li>Pipelines should ship with Nextflow profiles for singularity that pull from the docker repository.</li>
    </ul>
  </li>
  <li>Continuous integration testing</li>
  <li>Stable release tags</li>
  <li>Common pipeline structure and usage</li>
  <li>Excellent documentation</li>
  <li>A responsible contact person / GitHub username
    <ul>
      <li>This will typically be the main person behind the pipeline development</li>
      <li>This person should be responsible for basic maintenance and questions</li>
    </ul>
  </li>
  <li>The pipeline must not have any failures in the <code class="highlighter-rouge">nf-core lint</code> tests
    <ul>
      <li>These tests are run by the <a href="https://github.com/nf-core/tools">nf-core/tools</a> package and validate the requirements listed on this page.</li>
      <li>You can see the list of tests and how to pass them on the <a href="errors">error codes page</a>.</li>
    </ul>
  </li>
</ul>

<p><em>(any point can be skipped, given a good enough reason…)</em></p>

<p>If possible, it’s great if pipelines can also have:</p>

<ul>
  <li>Software bundled using <a href="https://bioconda.github.io/">bioconda</a>
    <ul>
      <li>See below for how to use a bioconda env script automatically build docker and singularity containers, meaning you have a single file to maintain.</li>
    </ul>
  </li>
  <li>Optimised output file formats
    <ul>
      <li>Pipelines should generate <code class="highlighter-rouge">CRAM</code> alignment files by default, but have a <code class="highlighter-rouge">--bam</code> option to generate <code class="highlighter-rouge">BAM</code> outputs if required by the user.</li>
    </ul>
  </li>
  <li>Explicit support for running in cloud environments
    <ul>
      <li>For example, use of <a href="https://ewels.github.io/AWS-iGenomes/">AWS-iGenomes</a></li>
    </ul>
  </li>
  <li>Benchmarks from running on cloud environments such as <a href="https://aws.amazon.com/">AWS</a></li>
</ul>

<h2 id="pipeline-organisation">Pipeline organisation</h2>
<p>It’s highly recommended that pipelines are built using the <a href="https://github.com/nf-core/cookiecutter">cookiecutter</a> starter template, as future developments are likely to be based on this assumption (see <a href="#plans-for-the-future"><em>future plans</em></a> below).</p>

<h2 id="coding-style">Coding style</h2>

<p>Pipelines must:</p>

<ul>
  <li>Use config profiles to organise anything that is not common for all users</li>
  <li>Run with as little input as possible
    <ul>
      <li>Metadata (eg. be able to run with just FastQ files, where possible)</li>
      <li>Reference files (eg. auto-generate missing reference files, where possible)</li>
    </ul>
  </li>
  <li>Keep only code for the latest stable on the main <code class="highlighter-rouge">master</code> branch.
    <ul>
      <li>The main development code should be kept in a branch called <code class="highlighter-rouge">development</code></li>
    </ul>
  </li>
  <li>Use GitHub releases and keep a detailed changelog file</li>
  <li><em>..and many more :)</em></li>
</ul>

<h1 id="how-to-add-a-new-pipeline">How to add a new pipeline</h1>

<ol>
  <li>Join the <a href="https://github.com/nf-core/nf-core.github.io/issues/1">nf-core GitHub organisation</a></li>
  <li>Create a pipeline repository in the organisation
    <ul>
      <li><em>If starting from scratch</em>
        <ul>
          <li>Ask an admin to create a new pipeline repository for you and add you as a collaborator.</li>
          <li>This is the best option, as this repository will then be recognised as the “head fork” by GitHub.</li>
        </ul>
      </li>
      <li><em>If you already have a pipeline</em>
        <ul>
          <li>Just fork your pipeline to the nf-core GitHub organisation</li>
        </ul>
      </li>
    </ul>
  </li>
  <li>Make sure that your pipeline <code class="highlighter-rouge">README.md</code> file has a big warning on the front saying that it’s under development</li>
  <li>When you’re happy with it, ping <a href="https://github.com/orgs/nf-core/teams/admin">@nf-core/admin</a> for a code review</li>
  <li>Once the pipeline has been approved, you can remove the development warning on the homepage and the pipeline will be added to the website.</li>
</ol>

<h1 id="plans-for-the-future">Plans for the future</h1>
<h2 id="base-change-automation">Base change automation</h2>
<p>In the future it would be great to have automated bots listen to changes in a base pipeline cookiecutter template. Changes (eg. in pull-requests) could be compared against all other pipelines; any that also apply elsewhere could then be flagged or even modified with automated pull-requests. The reverse could also apply - changes in a downstream pipeline that are shared in the core template could be acted upon.</p>

<h2 id="helper-scripts">Helper scripts</h2>
<p>It’s also possible that we could write a python package with helper scripts to make it easier to create custom config files, list available pipelines, check for updates, pull singularity images and so on. This could then be packaged using bioconda and the python package index as its own stand-alone tool.</p>

<p>Helper scripts could also be very useful when developing and testing pipelines. For example, linting required features and code style.</p>




<?php include('../includes/footer.php'); ?>
