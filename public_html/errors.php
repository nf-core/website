<?php
$title = 'Errors';
$subtitle = 'Find out details about different error codes from the <code>nf-core lint</code> command.';
include('../includes/header.php');
?>

<h1 id="linting-errors">Linting Errors</h1>

<p>This page contains detailed descriptions of the tests done by the <a href="https://github.com/nf-core/tools">nf-core/tools</a> package. Linting errors should show URLs next to any failures that link to the relevant heading below.</p>

<h2 id="1">Error #1 - File not found</h2>
<p>nf-core pipelines should adhere to a common file structure for consistency. The lint test looks for the following required files:</p>

<ul>
  <li><code>nextflow.config</code>
    <ul>
      <li>The main nextflow config file</li>
    </ul>
  </li>
  <li><code>Dockerfile</code>
    <ul>
      <li>A docker build script to generate a docker image with the required software</li>
    </ul>
  </li>
  <li><code>.travis.yml</code> or <code>circle.yml</code>
    <ul>
      <li>A config file for automated continuous testing with either <a href="https://travis-ci.org/">Travis CI</a> or <a href="https://circleci.com/">Circle CI</a></li>
    </ul>
  </li>
  <li><code>LICENSE</code>, <code>LICENSE.md</code>, <code>LICENCE.md</code> or <code>LICENCE.md</code>
    <ul>
      <li>The MIT licence. Copy from <a href="https://raw.githubusercontent.com/nf-core/tools/master/LICENSE">here</a>.</li>
    </ul>
  </li>
  <li><code>README.md</code>
    <ul>
      <li>A well written readme file in markdown format</li>
    </ul>
  </li>
  <li><code>CHANGELOG.md</code>
    <ul>
      <li>A markdown file listing the changes for each pipeline release</li>
    </ul>
  </li>
  <li><code>docs/README.md</code>, <code>docs/output.md</code> and <code>docs/usage.md</code>
    <ul>
      <li>A <code>docs</code> directory with an index <code>README.md</code>, usage and output documentation</li>
    </ul>
  </li>
</ul>

<p>The following files are suggested but not a hard requirement. If they are missing they trigger a warning:</p>

<ul>
  <li><code>main.nf</code>
    <ul>
      <li>It’s recommended that the main workflow script is called <code>main.nf</code></li>
    </ul>
  </li>
  <li><code>conf/base.config</code>
    <ul>
      <li>A <code>conf</code> directory with at least one config called <code>base.config</code></li>
    </ul>
  </li>
  <li><code>tests/run_test.sh</code>
    <ul>
      <li>A bash script to run the pipeline test run</li>
    </ul>
  </li>
</ul>

<h2 id="2">Error #2 - Docker file check failed</h2>
<p>Pipelines should have a file called <code>Dockerfile</code> in their root directory.
This is used for automated docker image builds. This test checks that the file
exists and contains at least the string <code>FROM </code>.</p>

<h2 id="3">Error #3 - Licence check failed</h2>
<p>nf-core pipelines must ship with an open source <a href="https://choosealicense.com/licenses/mit/">MIT licence</a>.</p>

<p>This test fails if the following conditions are not met:</p>

<ul>
  <li>No licence file found
    <ul>
      <li><code>LICENSE</code>, <code>LICENSE.md</code>, <code>LICENCE.md</code> or <code>LICENCE.md</code></li>
    </ul>
  </li>
  <li>Licence file contains fewer than 4 lines of text</li>
  <li>File does not contain the string <code>without restriction</code></li>
  <li>Licence contains template placeholders
    <ul>
      <li><code>[year]</code>, <code>[fullname]</code>, <code>&lt;YEAR&gt;</code>, <code>&lt;COPYRIGHT HOLDER&gt;</code>, <code>&lt;year&gt;</code> or <code>&lt;copyright holders&gt;</code></li>
    </ul>
  </li>
</ul>

<h2 id="4">Error #4 - Nextflow config check failed</h2>
<p>nf-core pipelines are required to be configured with a minimal set of variable
names. This test fails or throws warnings if required variables are not set.</p>

<blockquote>
  <p><strong>Note:</strong> These config variables must be set in <code>nextflow.config</code> or another config
file imported from there. Any variables set in nextflow script files (eg. <code>main.nf</code>)
are not checked and will be assumed to be missing.</p>
</blockquote>

<p>The following variables fail the test if missing:</p>

<ul>
  <li><code>params.version</code>
    <ul>
      <li>The version of this pipeline. This should correspond to a <a href="https://help.github.com/articles/creating-releases/">GitHub release</a>.</li>
    </ul>
  </li>
  <li><code>params.nf_required_version</code>
    <ul>
      <li>The minimum version of Nextflow required to run the pipeline.</li>
      <li>This should correspond to the <code>NXF_VER</code> version tested by Travis.</li>
    </ul>
  </li>
  <li><code>params.outdir</code>
    <ul>
      <li>A directory in which all pipeline results should be saved</li>
    </ul>
  </li>
  <li><code>manifest.description</code>
    <ul>
      <li>A description of the pipeline</li>
    </ul>
  </li>
  <li><code>manifest.homePage</code>
    <ul>
      <li>The homepage for the pipeline. Should be the nf-core GitHub repository URL,
so beginning with <code>https://github.com/nf-core/</code></li>
    </ul>
  </li>
  <li><code>timeline.enabled</code>, <code>trace.enabled</code>, <code>report.enabled</code>, <code>dag.enabled</code>
    <ul>
      <li>The nextflow timeline, trace, report and DAG should be enabled by default</li>
    </ul>
  </li>
  <li><code>process.cpus</code>, <code>process.memory</code>, <code>process.time</code>
    <ul>
      <li>Default CPUs, memory and time limits for tasks</li>
    </ul>
  </li>
</ul>

<p>The following variables throw warnings if missing:</p>

<ul>
  <li><code>manifest.mainScript</code>
    <ul>
      <li>The filename of the main pipeline script (recommended to be <code>main.nf</code>)</li>
    </ul>
  </li>
  <li><code>timeline.file</code>, <code>trace.file</code>, <code>report.file</code>, <code>dag.file</code>
    <ul>
      <li>Default filenames for the timeline, trace and report</li>
      <li>Should be set to a results folder, eg: <code>${params.outdir}/pipeline_info/trace.[workflowname].txt"</code>`</li>
      <li>The DAG file path should end with <code>.svg</code>
        <ul>
          <li>If Graphviz is not installed, Nextflow will generate a <code>.dot</code> file instead</li>
        </ul>
      </li>
    </ul>
  </li>
  <li><code>process.container</code>
    <ul>
      <li>A single default container for use by all processes</li>
    </ul>
  </li>
  <li><code>params.reads</code>
    <ul>
      <li>Input parameter to specify input data (typically FastQ files / pairs)</li>
    </ul>
  </li>
  <li><code>params.singleEnd</code>
    <ul>
      <li>Specify to work with single-end sequence data instead of default paired-end</li>
      <li>Used with Nextflow: <code>.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )</code></li>
    </ul>
  </li>
</ul>

<h2 id="5">Error #5 - Continuous Integration configuration</h2>
<p>nf-core pipelines must have CI testing with Travis or Circle CI.</p>

<p>This test fails if the following happens:</p>

<ul>
  <li><code>.travis.yml</code> does not contain the string <code>nf-core lint ${TRAVIS_BUILD_DIR}</code> under <code>script</code></li>
  <li><code>.travis.yml</code> does not test the Nextflow version specified in the pipeline as <code>nf_required_version</code>
    <ul>
      <li>This is expected in the <code>env</code> section of the config, eg:
        <pre><code>env:
    - NXF_VER=0.27.0
    - NXF_VER=’’</code></pre>
      </li>
      <li>At least one of these <code>NXF_VER</code> variables must match the <code>params.nf_required_version</code> version specified in the pipeline config</li>
      <li>Other variables can be specified on these lines as long as they are space separated.</li>
    </ul>
  </li>
</ul>

<h2 id="6">Error #6 - Repository <code>README.md</code> tests</h2>
<p>The <code>README.md</code> files for a project are very important and must meet some requirements:</p>

<ul>
  <li>Nextflow badge
    <ul>
      <li>If no Nextflow badge is found, a warning is given</li>
      <li>If a badge is found but the version doesn’t match the minimum version in the config file, the test fails</li>
      <li>Example badge code:
        <div class="language-markdown highlighter-rouge"><div class="highlight"><pre class="highlight"><code>  <span class="p">[</span><span class="nv">![Nextflow</span><span class="p">](</span><span class="sx">https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg</span><span class="p">)</span>](https://www.nextflow.io/)
</code></pre></div>        </div>
      </li>
    </ul>
  </li>
  <li>Bioconda badge
    <ul>
      <li>If your pipeline contains a file called <code>environment.yml</code>, a bioconda badge is required</li>
      <li>Required badge code:
        <div class="language-markdown highlighter-rouge"><div class="highlight"><pre class="highlight"><code>  <span class="p">[</span><span class="nv">![install with bioconda</span><span class="p">](</span><span class="sx">https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg</span><span class="p">)</span>](http://bioconda.github.io/)
</code></pre></div>        </div>
      </li>
    </ul>
  </li>
</ul>

<h2 id="7">Error #7 - Pipeline and container version numbers</h2>

<blockquote>
  <p>This test only runs when <code>--release</code> is set or <code>$TRAVIS_BRANCH</code> is equal to <code>master</code></p>
</blockquote>

<p>These tests look at <code>params.container</code>, <code>process.container</code> and <code>$TRAVIS_TAG</code>, only
if they are set.</p>

<ul>
  <li>Container name must have a tag specified (eg. <code>nfcore/pipeline:version</code>)</li>
  <li>Container tag / <code>$TRAVIS_TAG</code> must contain only numbers and dots</li>
  <li>Tags and <code>$TRAVIS_TAG</code> must all match one another</li>
</ul>

<h2 id="8">Error #8 - Conda environment tests</h2>

<blockquote>
  <p>These tests only run when your pipeline has a root file called <code>environment.yml</code></p>
</blockquote>

<ul>
  <li>The environment <code>name</code> must match the pipeline name and version
    <ul>
      <li>The pipeline name is found from the Nextflow config <code>manifest.homePage</code>,
which assumes that the URL is in the format <code>github.com/nf-core/[pipeline-name]</code></li>
      <li>Example: For <code>github.com/nf-core/test</code> version 1.4, the conda environment name should be <code>nfcore-test-1.4</code></li>
    </ul>
  </li>
</ul>

<p>Each dependency is checked using the <a href="https://api.anaconda.org/docs">Anaconda API service</a>.
Dependency sublists are ignored with the exception of <code>- pip</code>: these packages are also checked
for pinned version numbers and checked using the <a href="https://wiki.python.org/moin/PyPIJSON">PyPI JSON API</a>.</p>

<p>Note that conda dependencies with pinned channels (eg. <code>conda-forge::openjdk</code>) are fine
and should be handled by the linting properly.</p>

<p>Each dependency can have the following lint failures and warnings:</p>

<ul>
  <li>(Test failure) Dependency does not have a pinned version number, eg. <code>toolname=1.6.8</code></li>
  <li>(Test failure) The package cannot be found on any of the listed conda channels (or PyPI if <code>pip</code>)</li>
  <li>(Test warning) A newer version of the package is available</li>
</ul>

<h2 id="9">Error #9 - Dockerfile for use with Conda environments</h2>

<blockquote>
  <p>This test only runs if there is both <code>environment.yml</code>
and <code>Dockerfile</code> present in the workflow.</p>
</blockquote>

<p>If a workflow has a conda <code>environment.yml</code> file (see above), the <code>Dockerfile</code> should use this
to create the container. Such <code>Dockerfile</code>s can usually be very short, eg:</p>

<pre><code>FROM nfcore/base
LABEL authors="your@email.com" \
      description="Docker image containing all requirements for nf-core/EXAMPLE pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml &amp;&amp; conda clean -a
ENV PATH /opt/conda/envs/nfcore-EXAMPLE/bin:$PATH
</code></pre>

<p>To enforce this minimal <code>Dockerfile</code> and check for common copy+paste errors, we require
that the above template is used.
Failures are generated if the <code>FROM</code>, <code>COPY</code>, <code>RUN</code> and <code>ENV</code> statements above are not present.
These lines must be an exact copy of the above example, with the exception that
the <code>ENV PATH</code> must reference the name of your pipeline instead of <code>nfcore-EXAMPLE</code>.</p>

<p>Additional lines and different metadata can be added without causing the test to fail.</p>



<?php include('../includes/footer.php'); ?>
