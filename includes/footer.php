<?php if(isset($title) and $title): ?>
      </div>
    </div>
 <?php endif; ?>


    <footer>
      <div class="container">
        <div class="row">
          <div class="col-6 col-md-4">
            <a class="text-muted" href="/"><img height="30px" src="assets/img/logo/nf-core-logo.svg" /></a>
            <small class="d-block mb-3 text-muted">Making awesome workflows since &copy; 2018</small>
            <small class="d-block mb-3 text-muted">
              Website by <a href="http://phil.ewels.co.uk/" class="text-muted">Phil Ewels</a>.
              Icons from <a href="http://www.flaticon.com/" class="text-muted">flaticon.com</a>
              and <a href="https://fontawesome.com/" class="text-muted">fontawesome.com</a>.
              <a href="http://getbootstrap.com/" class="text-muted">Bootstrap</a> CSS framework,
              <a href="http://jquery.com/" class="text-muted">jQuery</a> JS and syntax colouring
              with <a href="https://highlightjs.org/" class="text-muted">highlight.js</a>.
            </small>
          </div>
          <div class="col-6 col-md">
            <h5>Getting Started</h5>
            <ul class="list-unstyled">
              <li><a class="text-muted" href="#">Workflow features</a></li>
              <li><a class="text-muted" href="#">Using nextflow</a></li>
              <li><a class="text-muted" href="/pipelines">Available pipelines</a></li>
            </ul>
          </div>
          <div class="col-6 col-md">
            <h5>For Authors</h5>
            <ul class="list-unstyled">
              <li><a class="text-muted" href="/developer_docs">Documentation</a></li>
              <li><a class="text-muted" href="/errors">Error codes</a></li>
              <li><a class="text-muted" href="/tools">Helper tools</a></li>
            </ul>
          </div>
          <div class="col-6 col-md">
            <h5>About nf-core</h5>
            <ul class="list-unstyled">
              <li><a class="text-muted" href="/about#history">Project history</a></li>
              <li><a class="text-muted" href="/about#contributors">List of contributors</a></li>
              <li><a class="text-muted" href="/developer_docs#getting-help">Where to get help</a></li>
            </ul>
          </div>
        </div>
      </div>
    </footer>


    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
    <script>window.jQuery || document.write('assets/js/jquery-slim.min.js"><\/script>')</script>
    <script src="assets/js/popper.min.js"></script>
    <script src="assets/js/bootstrap.min.js"></script>
    <script src="assets/js/highlight.pack.js"></script>
    <script src="assets/js/nf-core.js?c=<?php echo $git_sha; ?>"></script>
  </body>
</html>
