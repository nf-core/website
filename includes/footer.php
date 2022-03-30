<?php
if (isset($title) and $title) {
  echo '</div></div>';
}

// Subfooter, if given
if (isset($md_github_url) and $md_github_url) {
  $subfooter = '<p class="mb-0"><i class="fab fa-github"></i> Read this page on GitHub: <code><a href="' . $md_github_url . '">' . $md_github_url . '</a></code>';
}
if (isset($subfooter) and $subfooter) {
  echo '<footer class="subfooter"><div class="container">' . $subfooter . '</div></footer>';
}
?>

<footer class="footer">
  <div class="container">
    <div class="row">
      <div class="col-lg-3 mb-3">
        <a href="/">
          <img height="30px" src="/assets/img/logo/nf-core-logo.svg" class="hide-dark">
          <img height="30px" src="/assets/img/logo/nf-core-logo-darkbg.svg" class="hide-light">
        </a>
        <small class="d-block mb-3">Making awesome workflows since 2018</small>
        <small class="d-flex mb-3">Supported by<div class="social-icons">
            <a href="/about#czi-eoss">
              <img src="/assets/img/contributors-colour/CZI-alt.svg" alt="CZI" class="" style="max-width: 50px">
            </a>
            <span class="ms-1">and</span>
            <a href="/about#aws">
              <img src="/assets/img/contributors-white/aws.svg" alt="Amazon Web Services" class="hide-light" style="max-width: 50px">
              <img src="/assets/img/contributors-colour/aws.svg" alt="Amazon Web Services" class="hide-dark" style="max-width: 50px">
            </a>
          </div>
        </small>
        <small class="d-block mb-3">
          See the source code for this website on GitHub: <a href="https://github.com/nf-core/nf-co.re" target="_blank">https://github.com/nf-core/nf-co.re</a>
        </small>
        <div class="d-md-flex d-print-none">
          <div class="btn-toolbar mb-3 me-4" role="toolbar">
            <div class="theme-switcher btn-group btn-group-sm" role="group">

              <input type="radio" class="btn-check" id="theme-auto" name="theme-auto" value="auto" autocomplete="off" <?php if ($theme == 'auto') echo 'checked'; ?>>
              <label class="btn btn-secondary" for="theme-auto" data-bs-toggle=" tooltip" title="Auto Light / Dark"><i class="fas fa-adjust"></i></label>

              <input type="radio" class="btn-check" id="theme-light" name="theme-light" value="light" autocomplete="off" <?php if ($theme == 'light') echo 'checked'; ?>>
              <label class="btn btn-secondary" for="theme-light" data-bs-toggle=" tooltip" title="Light Theme"><i class="fas fa-sun"></i></label>

              <input type="radio" class="btn-check" id="theme-dark" name="theme-dark" value="dark" autocomplete="off" <?php if ($theme == 'dark') echo 'checked'; ?>>
              <label class="btn btn-secondary" for="theme-dark" data-bs-toggle=" tooltip" title="Dark Theme"><i class="fas fa-moon"></i></label>
            </div>
          </div>

          <div class="social-icons mb-3 d-print-none">
            <a href="https://nfcore.slack.com/" target="_blank" title="Slack" data-bs-toggle="tooltip">
              <img src="/assets/img/slack.svg" />
            </a>
            <a href="https://github.com/nf-core/" target="_blank" class="social-github" title="GitHub" data-bs-toggle="tooltip">
              <img src="/assets/img/github.svg" />
            </a>
            <a href="https://twitter.com/nf_core" target="_blank" title="Twitter" data-bs-toggle="tooltip">
              <img src="/assets/img/twitter.svg" />
            </a>
            <a href="https://www.youtube.com/c/nf-core" target="_blank" title="YouTube" data-bs-toggle="tooltip">
              <img src="/assets/img/youtube.svg" />
            </a>
          </div>
        </div>

      </div>
      <div class="col-sm-6 col-lg-3 offset-lg-1 mb-3 d-print-none">
        <h5>Getting Started</h5>
        <ul class="list-unstyled">
          <li><a href="/pipelines">Available pipelines</a></li>
          <li><a href="/modules">Available modules</a></li>
          <li><a href="/tools">Helper tools</a></li>
          <li><a href="/usage/introduction">Getting started</a></li>
          <li><a href="/usage/installation">Installation</a></li>
          <li><a href="/usage/configuration">Pipeline configuration</a></li>
          <li><a href="/usage/offline">Running offline</a></li>
          <li><a href="/usage/usage_tutorials">Usage tutorials</a></li>
          <li><a href="/usage/reference_genomes">Reference genomes</a></li>
          <li><a href="/usage/data_management">Data Management</a></li>
          <li><a href="/usage/troubleshooting">Troubleshooting</a></li>
          <li><a href="/usage/nextflow">Nextflow resources</a></li>
        </ul>
      </div>
      <div class="col-sm-6 col-lg-3 mb-3 d-print-none">
        <h5>For Authors</h5>
        <ul class="list-unstyled">
          <li><a href="/developers/guidelines">Guidelines</a></li>
          <li><a href="/developers/adding_pipelines">Adding a new pipeline</a></li>
          <li><a href="/developers/modules">Adding a new module</a></li>
          <li><a href="/developers/release_checklist">Release checklist</a></li>
          <li><a href="/tools-docs">Lint error codes</a></li>
          <li><a href="/developers/sync">Template synchronisation</a></li>
          <li><a href="/developers/developer_tutorials">Developer tutorials</a></li>
          <li><a href="/developers/editor_plugins">Code editor plugins</a></li>
          <li><a href="/developers/design_guidelines">Graphic design guidelines</a></li>
        </ul>
      </div>
      <div class="col-sm-6 col-lg-2 mb-3 d-print-none">
        <h5>About nf-core</h5>
        <ul class="list-unstyled">
          <li><a href="/about">About nf-core</a></li>
          <li><a href="/events">Events</a></li>
          <li><a href="/community">Community</a></li>
          <li><a href="/stats">Statistics</a></li>
          <li><a href="/publications">Publications</a></li>
          <li><a href="/code_of_conduct">Code of conduct</a></li>
          <li><a href="/join" class="mt-1 btn btn-outline-success">Join nf-core</a></li>
        </ul>
      </div>
    </div>
  </div>
</footer>
<?php
if (isset($end_of_html) and $end_of_html) {
  echo $end_of_html;
}
?>
</body>

</html>