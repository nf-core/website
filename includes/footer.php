<?php
if(isset($title) and $title){
  echo '</div></div>';
}

// Subfooter, if given
if(isset($md_github_url) and $md_github_url){
  $subfooter = '<p class="mb-0"><i class="fab fa-github"></i> Read this page on GitHub: <code><a href="'.$md_github_url.'">'.$md_github_url.'</a></code>';
}
if(isset($subfooter) and $subfooter){
  echo '<footer class="subfooter"><div class="container">'.$subfooter.'</div></footer>';
}
?>

    <footer class="footer">
      <div class="container">
        <div class="row">
          <div class="col-6 col-md-4">
            <a href="/">
                <img height="30px" src="/assets/img/logo/nf-core-logo.svg" class="hide-dark">
                <img height="30px" src="/assets/img/logo/nf-core-logo-darkbg.svg" class="hide-light">
            </a>
            <small class="d-block mb-3">Making awesome workflows since 2018</small>
            <small class="d-block mb-3">
              Website by <a href="http://phil.ewels.co.uk/">Phil Ewels</a>.
              Icons from <a href="http://www.flaticon.com/">flaticon.com</a>
              and <a href="https://fontawesome.com/">fontawesome.com</a>.
              <a href="http://getbootstrap.com/">Bootstrap</a> CSS framework,
              <a href="http://jquery.com/">jQuery</a> JS and syntax colouring
              with <a href="https://highlightjs.org/">highlight.js</a>.
            </small>
            <div class="d-md-flex">
              <div class="btn-toolbar mb-3 mr-4" role="toolbar">
                <div class="theme-switcher border btn-group btn-group-sm btn-group-toggle" data-toggle="buttons">
                  <label class="btn btn-light <?php if($theme == 'auto') echo 'active'; ?>" data-toggle="tooltip" title="Auto Light / Dark">
                    <input type="radio" value="auto" autocomplete="off" checked> <i class="fas fa-adjust"></i>
                  </label>
                  <label class="btn btn-light <?php if($theme == 'light') echo 'active'; ?>" data-toggle="tooltip" title="Light Theme">
                    <input type="radio" value="light" autocomplete="off"> <i class="fas fa-sun"></i>
                  </label>
                  <label class="btn btn-light <?php if($theme == 'dark') echo 'active'; ?>" data-toggle="tooltip" title="Dark Theme">
                    <input type="radio" value="dark" autocomplete="off"> <i class="fas fa-moon"></i>
                  </label>
                </div>
              </div>

              <div class="social-icons mb-3">
                <a href="https://nfcore.slack.com/" target="_blank" title="Slack" data-toggle="tooltip">
                 <img src="/assets/img/slack.svg" />
                </a>
                <a href="https://github.com/nf-core/" target="_blank" class="social-github" title="GitHub" data-toggle="tooltip">
                 <img src="/assets/img/github.svg" />
                </a>
                <a href="https://twitter.com/nf_core" target="_blank" title="twitter" data-toggle="tooltip">
                 <img src="/assets/img/twitter.svg" />
                </a>
                <a href="https://www.youtube.com/c/nf-core" target="_blank" title="YouTube" data-toggle="tooltip">
                 <img src="/assets/img/youtube.svg" />
                </a>
              </div>
            </div>

          </div>
          <div class="col-6 col-md">
            <h5>Getting Started</h5>
            <ul class="list-unstyled">
              <li><a href="/pipelines">Available pipelines</a></li>
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
          <div class="col-6 col-md">
            <h5>For Authors</h5>
            <ul class="list-unstyled">
              <li><a href="/developers/guidelines">Guidelines</a></li>
              <li><a href="/developers/adding_pipelines">Adding a new pipeline</a></li>
              <li><a href="/developers/release_checklist">Release checklist</a></li>
              <li><a href="/tools-docs">Lint error codes</a></li>
              <li><a href="/developers/sync">Template synchronisation</a></li>
              <li><a href="/developers/developer_tutorials">Developer tutorials</a></li>
              <li><a href="/developers/editor_plugins">Code editor plugins</a></li>
              <li><a href="/developers/design_guidelines">Graphic design guidelines</a></li>
            </ul>
          </div>
          <div class="col-6 col-md">
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
    if(isset($end_of_html) and $end_of_html){
      echo $end_of_html;
    }
    ?>
  </body>
</html>
