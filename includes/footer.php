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
            <a href="/"><img height="30px" src="/assets/img/logo/nf-core-logo.svg" /></a>
            <small class="d-block mb-3">Making awesome workflows since &copy; 2018</small>
            <small class="d-block mb-3">
              Website by <a href="http://phil.ewels.co.uk/">Phil Ewels</a>.
              Icons from <a href="http://www.flaticon.com/">flaticon.com</a>
              and <a href="https://fontawesome.com/">fontawesome.com</a>.
              <a href="http://getbootstrap.com/">Bootstrap</a> CSS framework,
              <a href="http://jquery.com/">jQuery</a> JS and syntax colouring
              with <a href="https://highlightjs.org/">highlight.js</a>.
            </small>
          </div>
          <div class="col-6 col-md">
            <h5>Getting Started</h5>
            <ul class="list-unstyled">
              <li><a href="usage/introduction">Using nf-core</a></li>
              <li><a href="/pipelines">Available pipelines</a></li>
              <li><a href="/tools">Helper tools</a></li>
              <li><a href="usage/nextflow_tutorial">Nextflow tutorial</a></li>
            </ul>
          </div>
          <div class="col-6 col-md">
            <h5>For Authors</h5>
            <ul class="list-unstyled">
              <li><a href="developers/guidelines">Guidelines</a>
              <li><a href="developers/adding_pipelines">Adding a new pipeline</a>
              <li><a href="/errors">Lint error codes</a>
              <li><a href="developers/sync">Template synchronisation</a>
            </ul>
          </div>
          <div class="col-6 col-md">
            <h5>About nf-core</h5>
            <ul class="list-unstyled">
              <li><a href="/about#contributors">List of contributors</a></li>
              <li><a href="/about#history">Project history</a></li>
              <li><a href="/about#contact">Get in touch</a></li>
            </ul>
          </div>
        </div>
      </div>
    </footer>

  </body>
</html>
