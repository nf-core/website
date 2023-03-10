<?php

require '../vendor/autoload.php'; // for spyc

// Get a random subset of contributor institute logos
$contributors = spyc_load_file('../nf-core-contributors.yaml');
$contributors_img_list = [];
foreach ($contributors['contributors'] as $idx => $c) {
    if (isset($c['image_fn']) and $c['image_fn'] and isset($c['full_name']) and $c['full_name']) {
        $card_id = preg_replace('/[^a-z]+/', '-', strtolower($c['full_name']));
        $img_path = 'assets/img/contributors-white/' . $c['image_fn'];
        if (file_exists($img_path)) {
            $contributors_img_list[] =
                '<a href="/community#' .
                $card_id .
                '"><img src="' .
                $img_path .
                '" data-bs-placement="bottom" data-bs-toggle="tooltip" title="' .
                $c['full_name'] .
                '"></a>';
        }
    }
}
// Shuffle and truncate the list
shuffle($contributors_img_list);

//Check if there is an event
$md_github_url = 'https://github.com/sanger-tol/nf-co.re/tree/main/markdown/events';
$header_btn_url = 'https://nf-co.re/events/rss';

# To get parse_md_front_matter() and sanitise_date_meta() functions
require_once '../includes/functions.php';

if ($curr_event) {
    // Shared function to prep nicely formatted output
    $curr_event['meta'] = prep_current_event($curr_event);
    // Dropdown button to visit event
    $curr_event['meta']['location_dropdown'] = '';
    if (array_key_exists('location_url', $curr_event) && $curr_event['ongoing']) {
        if (count($curr_event['location_url']) == 1) {
            $url = $curr_event['location_url'][0];
            if ($url[0] == '#') {
                $url = $curr_event['url'] . $url;
            }
            $m = $curr_event['meta']['location_url_meta'];

            $curr_event['meta']['location_dropdown'] =
                '<a class="btn btn-success me-2 mb-2" href="' . $url . '">' . $m[0]['icon'] . ' Watch now</a>';
        } else {
            $curr_event['meta']['location_dropdown'] = '
        <div class="dropup me-2 mb-2">
          <a class="btn btn-success dropdown-toggle" href="#" role="button" data-bs-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
            Watch now
          </a>

          <div class="dropdown-menu">
          ';
            foreach ($curr_event['location_url'] as $idx => $url) {
                $m = $curr_event['meta']['location_url_meta'][$idx];
                if ($url[0] == '#') {
                    $url = $curr_event['url'] . $url;
                }
                $curr_event['meta']['location_dropdown'] .=
                    '<a class="dropdown-item" href="' .
                    $url .
                    '" target="_blank">' .
                    $m['icon'] .
                    ' <code>' .
                    $url .
                    '</code></a>' .
                    "\n";
            }
            $curr_event['meta']['location_dropdown'] .= '</div></div>';
        }
    }
    // Countdown timer for upcoming events
    if (!$curr_event['ongoing']) {
        $dtF = new \DateTime('@0');
        $dtT = new \DateTime('@' . (time() - $curr_event['start_ts']));
        $dtDiff = $dtF->diff($dtT);
        if ($dtDiff->format('%d') == '0') {
            $countdown_text = $dtDiff->format('%hh %Im %Ss');
        } else {
            $countdown_text = $dtDiff->format('%d days,<br>%hh %Im %Ss');
        }
        $curr_event['meta']['countdown'] =
            "
    <script type=\"text/javascript\">
        setInterval(function(){
          var eventTime = " .
            $curr_event['start_ts'] .
            " * 1000;
          var currentTime = Date.now();
          var delta = Math.abs(eventTime - currentTime) / 1000;

          var days = Math.floor(delta / 86400);
          delta -= days * 86400;

          var hours = Math.floor(delta / 3600) % 24;
          delta -= hours * 3600;

          var minutes = Math.floor(delta / 60) % 60;
          delta -= minutes * 60;

          var seconds = Math.floor(delta) % 60;

          $('.countdown').html((days > 0 ? days + ' day' + (days > 1 ? 's':'') + ',<br>' : '') + hours  + 'h ' + minutes + 'm ' + seconds +'s')
        }, 1000);
    </script>
    <h5 class=\"pt-4\">Event countdown:</h5>
    <p class=\"display-5 text-nowrap countdown\">" .
            $countdown_text .
            "</p>
    ";
    }
}

$import_moment = true;
include '../includes/header.php';
?>

<div class="homepage-header">
  <?php if ($curr_event): ?>
    <div class="mainpage-subheader-heading homepage-header-contents triangle-down pb-4">
      <div class="container-fluid text-start">
        <div class="row">
          <?php if ($curr_event['ongoing']): ?>
            <div class="col-lg-3 overflow-hidden d-none d-lg-block">
              <i class="fad fa-broadcast-tower homepage-header-fa-background"></i>
              <h4 class="display-4 pt-2">Ongoing event</h4>
            </div>
            <div class="col-lg-3 pt-2 pb-1 mb-2 overflow-hidden d-lg-none mainpage-subheader-heading-header">
              <h5 class="pt-2 font-weight-light">Ongoing event</h5>
            </div>
          <?php else: ?>
            <div class="col-lg-3 overflow-hidden d-none d-lg-block">
              <i class="fad fa-alarm-clock homepage-header-fa-background"></i>
              <h4 class="display-4 pt-2">Upcoming event</h4>
            </div>
            <div class="col-lg-3 pt-2 pb-1 mb-2 overflow-hidden d-lg-none mainpage-subheader-heading-header">
              <h5 class="pt-2 font-weight-light">Upcoming event</h5>
            </div>
          <?php endif; ?>
          <div class="col pt-lg-3 pb-lg-3 text-center text-lg-start">
            <h5 class="pt-2 pb-0 pb-lg-1"><a href="<?php echo $curr_event[
                'url'
            ]; ?>" class="text-success text-decoration-none"><?php echo $curr_event['title']; ?></a></h5>
            <p class="lead d-none d-sm-block"><a href="<?php echo $curr_event[
                'url'
            ]; ?>" class="text-body text-decoration-none"><?php echo $curr_event['subtitle']; ?></a></p>
            <p class="d-sm-none mb-2"><a href="<?php echo $curr_event[
                'url'
            ]; ?>" class="text-body text-decoration-none"><?php echo $curr_event['subtitle']; ?></a></p>
            <p class="d-none d-lg-block"><a href="<?php echo $curr_event[
                'url'
            ]; ?>" class="text-secondary text-decoration-none" <?php echo $curr_event['meta'][
    'nice_date_string'
][0]; ?>><?php echo $curr_event['meta'][
    'nice_date_string'
][1]; ?></a><span class="d-none d-lg-inline"> &nbsp; <?php echo $curr_event['meta']['event_type_badge']; ?></span></p>
            <p class="d-lg-none small mb-2"><a href="<?php echo $curr_event[
                'url'
            ]; ?>" class="text-secondary text-decoration-none" <?php echo $curr_event['meta'][
    'nice_date_string'
][0]; ?>><?php echo $curr_event['meta'][
    'nice_date_string'
][1]; ?></a><span class="d-none d-lg-inline"> &nbsp; <?php echo $curr_event['meta']['event_type_badge']; ?></span></p>
            <?php
            if ($curr_event['ongoing'] && isset($curr_event['youtube_embed'])): ?>
              <div class="btn-toolbar justify-content-center justify-content-lg-start">
                <?php echo $curr_event['meta']['location_dropdown']; ?>
                <a href="<?php echo $curr_event[
                    'url'
                ]; ?>" class="btn btn-outline-success mb-2 d-none d-lg-inline-block">Event Details</a>
                <a href="<?php echo $curr_event[
                    'url'
                ]; ?>" class="btn btn-sm btn-outline-success mb-2 d-lg-none">Event Details</a>
              </div>
            <?php endif;
            if (!$curr_event['ongoing']): ?>
              <a href="<?php echo $curr_event[
                  'url'
              ]; ?>" class="btn btn-outline-success mb-2 d-none d-lg-inline-block">Event Details</a>
              <a href="<?php echo $curr_event[
                  'url'
              ]; ?>" class="btn btn-sm btn-outline-success mb-2 d-lg-none">Event Details</a>
            <?php endif;
            ?>
          </div>
          <div class="col-lg-4 col-xl-3">
            <?php if ($curr_event['ongoing'] && isset($curr_event['youtube_embed'])):

                if (!is_array($curr_event['youtube_embed'])) {
                    $curr_event['youtube_embed'] = [$curr_event['youtube_embed']];
                }
                $video_id = get_youtube_id($curr_event['youtube_embed'][0]);
                ?>
              <div class="ratio ratio-16x9 mt-3 mb-5 d-none d-lg-block">
                <iframe width="560" height="315" src="https://www.youtube.com/embed/<?php echo $video_id; ?>" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
              </div>
              <?php if (count($curr_event['youtube_embed']) > 1): ?>
                <p class=""><a href="<?php echo $curr_event['url']; ?>" class="btn btn-success my-2">Watch <?php
echo count($curr_event['youtube_embed']) - 1;
d;
endif;
            elseif ($curr_event['ongoing']): ?>
              <div class="pt-lg-5">
                <?php echo $curr_event['meta']['location_dropdown']; ?>
                <a href="<?php echo $curr_event['url']; ?>" class="btn btn-outline-success mb-2">Event Details</a>
              </div>
            <?php else: ?>
              <div class="d-none d-lg-block">
                <?php echo $curr_event['meta']['countdown']; ?>
              </div>
            <?php endif; ?>
          </div>
        </div>
      </div>
    </div>
  <?php endif; ?>
  <div class="homepage-header-contents col-md-5">
    <h1>
      <img src="assets/img/logo/sanger-tol-logo.svg" class="hide-dark">
      <img src="assets/img/logo/sanger-tol-logo-darkbg.svg" class="hide-light">
    </h1>
    <p class="lead font-weight-normal">Workflows and tools to investigate the diversity of complex organisms.</p>
    <div class="homepage-cta-flex mb-5">
      <a class="homepage-cta" href="/pipelines">View Pipelines</a>
    </div>
  </div>
</div>
<div class="d-flex justify-content-center">
  <form class="searchbar_form homepage_search" action="search" method="get">
    <div class="input-group">
      <input type="search" class="form-control" placeholder="Search" name="q" required>
      <button type="submit" class="btn btn-outline-success">Search</button>
    </div>
  </form>
</div>

<div class="homepage-intro hexagon">
  <div class="container">
    <div class="row justify-content-center">
      <div class="col-md-4 col-xl-3 py-3 py-md-5 px-md-4">
        <i class="fad fa-industry-alt fa-3x mb-2"></i>
        <h3 class="text-white">For facilities</h3>
        <p class="lead">Highly optimised pipelines with excellent reporting. Validated releases ensure reproducibility.</p>
      </div>
      <div class="col-md-4 col-xl-3 py-3 py-md-5 px-md-4">
        <i class="fad fa-users fa-3x mb-2"></i>
        <h3 class="text-white">For users</h3>
        <p class="lead">Portable, documented and easy&nbsp;to&nbsp;use workflows.<br>Pipelines that you can trust.</p>
      </div>
      <div class="col-md-4 col-xl-3 py-3 py-md-5 px-md-4">
        <i class="fad fa-laptop-code fa-3x mb-2"></i>
        <h3 class="text-white">For developers</h3>
        <p class="lead">Companion templates and tools help to validate your code and simplify common tasks.</p>
      </div>
    </div>
    <div>
      <p class="d-flex position-relative justify-content-center text-white m-0 pb-3">
        <img src="/assets/img/pnas-logo.png" alt="PNAS" width="60px" class="float-start me-3">
        <a class="text-white stretched-link text-decoration-none d-flex flex-column justify-content-center" href="https://doi.org/10.1073/pnas.2115642118" target="_blank">
          Sequence locally, think globally: The Darwin Tree of Life Project
          <span class="underline m-auto"><em>PNAS</em> <strong>119 (4)</strong> e2115642118 (2022).
            <img src="/assets/img/OpenAccess.svg" alt="PNAS" height="20px" class="ms-2"></span>
        </a>
      </p>
    </div>
  </div>
</div>


<div class="container py-5 text-center lead" id="features">
  <p>Tree of Life production suite is designed to handle nature's diversity.</p>
  <p><strong>sanger-tol</strong> pipelines use Nextflow DSL2 and adhere to strict nf-core guidelines.</p>
</div>

<div id="features" class="container homepage-feature-boxes pb-5">
  <h3 class="mb-4 text-center">Fully featured pipelines</h3>
  <div class="row row-cols-1 row-cols-md-3 g-4">
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">Documentation</h5>
            <p class="card-text">Extensive documentation covering installation, usage and description of output files
              ensures that you won't be left in the dark.</p>
          </div>
          <img class="ms-3" height="100px" src="assets/img/docs.svg" />
        </div>
      </div>
    </div>
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">CI Testing</h5>
            <p class="card-text">Every time a change is made to the pipeline code,
              sanger-tol pipelines use continuous-integration testing to ensure that nothing has broken.</p>
          </div>
          <img class="ms-3" height="90px" src="assets/img/github-actions.svg" />
        </div>
      </div>
    </div>
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">Stable Releases</h5>
            <p class="card-text">sanger-tol pipelines use GitHub releases to tag stable versions of the code
              and software, making pipeline runs totally reproducible.</p>
          </div>
          <img class="ms-3" height="100px" src="assets/img/releases.svg" />
        </div>
      </div>
    </div>
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">Packaged software</h5>
            <p class="card-text">Pipeline dependencies are automatically downloaded and
              handled using Docker, Singularity, Conda or others.
              No need for any software installations.</p>
          </div>
          <i class="fad fa-box-alt fa-5x text-info ms-3"></i>
        </div>
      </div>
    </div>
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">Portable and reproducible</h5>
            <p class="card-text">
              Pipelines follow best-practices to ensure maximum
              portability and reproducibility.
              Scale and diversity in Tree of Life makes the pipelines 
              exceptionally well tested and easy to run.
            </p>
          </div>
          <i class="fas fa-server fa-5x text-success ms-3"></i>
        </div>
      </div>
    </div>
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body d-flex align-items-center">
          <div>
            <h5 class="card-title">Robust</h5>
            <p class="card-text">Pipelines are tested on data from different parts of the tree of life.
              sanger-tol pipelines are for all organisms not just model organisms.
            </p>
          </div>
          <img class="ms-3 hide-light" width="100px" src="/assets/img/contributors-white/aws.svg" alt="Amazon Web Services">
          <img class="ms-3 hide-dark" width="100px" src="/assets/img/contributors-colour/aws.svg" alt="Amazon Web Services">
        </div>
      </div>
    </div>
  </div>
</div>


<div class="bg-secondary pt-4 pb-5" style="--bs-bg-opacity: .1;">
  <div id="developers" class="container homepage-feature-boxes">
    <h3 class="mb-4 text-center">Developers: Not just another registry</h3>
    <div class="row row-cols-1 row-cols-md-3 g-4">
      <div class="cold">
        <div class="card h-100 shadow-sm">
          <div class="card-body">
            <h5 class="card-title">Develop <u>with</u> the community</h5>
            <div class='d-flex align-items-center'>
              <p class="card-text">Come and talk to us <em>before</em> you start writing a pipeline
                to find collaborators and check that your pipeline is suitable for sanger-tol.</p>
              <i class="fad fa-people-carry fa-5x text-secondary ms-3"></i>
            </div>
            <a href="join" class="btn btn-sm btn-outline-success arrow-hover"><span>Join sanger-tol</span></a>
          </div>
        </div>
      </div>
      <div class="cold">
        <div class="card h-100 shadow-sm">
          <div class="card-body">
            <h5 class="card-title">Start from the template</h5>
            <div class='d-flex align-items-center'>
              <p class="card-text">All pipelines and modules must be based on nf-core templates and be created using <a href="https://github.com/nf-core/tools">nf-core tools</a>.</p>
              </p>
              <i class="fas fa-magic fa-5x text-secondary ms-3"></i>
            </div>
            <a href="contributing/adding_pipelines" class="btn btn-sm btn-outline-success arrow-hover"><span>Read the docs</span></a>
          </div>
        </div>
      </div>
      <div class="cold">
        <div class="card h-100 shadow-sm">
          <div class="card-body">
            <h5 class="card-title">Collaborate, don't duplicate</h5>
            <div class='d-flex align-items-center'>
              <p class="card-text">We only allow one pipeline per analysis type.
                If a similar pipeline exists we'll ask you to add to that instead of making a new workflow.</p>
              <i class="fad fa-code-merge fa-5x text-secondary ms-3"></i>
            </div>
            <a href="contributing/guidelines" class="btn btn-sm btn-outline-success arrow-hover"><span>See the guidelines</span></a>
          </div>
        </div>
      </div>
    </div>
  </div>
</div>

<div class="bg-secondary py-5">
  <div class="container">

    <div class="row videos-row">
      <div class="col-lg-6">
        <div class="ratio ratio-16x9 hidden-xs hidden-sm">
          <iframe id="dtol-intro-video"  src="https://www.youtube.com/embed/aK1Ek39z4sA" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>
      </div>
      <div class="col-lg-6">
        <ul class="list-switch left video-chooser">
          <li><a href="https://youtu.be/aK1Ek39z4sA" data-src="https://www.youtube.com/embed/aK1Ek39z4sA" class="active"> Introduction to Darwin Tree of Life <em>(0:44)</em></a></li>
          <li><a href="https://youtu.be/Zort8tv7iRI" data-src="https://www.youtube.com/embed/Zort8tv7iRI"> Nextflow Summit: sanger-tol Production Engine <em>(14:57)</em></a></li>
          <li><a href="https://youtu.be/gUM9acK25tQ" data-src="https://www.youtube.com/embed/gUM9acK25tQ"> Introduction to nf-core <em>(1:01)</em></a></li>
        </ul>
      </div>
    </div>

  </div>
</div>

<div class="bg-dark py-5">
  <div class="container">
    <div class="row">
      <div class="col-sm-6">
        <h2 id="get-started" class="text-white">Get started in minutes</h2>
        <p class="text-white-50">Nextflow lets you run sanger-tol pipelines on virtually any computing environment.</p>
        <p class="text-white-50">All sanger-tol genomics pipelines come with built-in support for
          <a href="https://github.com/nf-core/tools" target="_blank" style="white-space:nowrap;">nf-core tools</a>.</p>
      </div>
      <div class="col-sm-6">
        <div class="card text-white bg-secondary border-secondary" style="--bs-bg-opacity: .3;">
          <div class="card-body">
            <pre class="text-white mb-0">
<span class="text-white-50"># Install nextflow</span>
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

<span class="text-white-50"># Launch the readmapping pipeline</span>
nextflow run sanger-tol/readmapping \
    --input samplesheet.csv \
    --fasta genome.fasta \
    -profile singularity

<span class="text-white-50"># Install nf-core tools</span>
pip install nf-core
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
      <a class="btn btn-lg btn-success arrow-hover" href="/pipelines"><span><i class="fas fa-toolbox me-1"></i> Available pipelines</span></a>
    </div>
    <div class="col-md-4 px-md-5 mb-5 mb-md-0">
      <h3>Run a pipeline</h3>
      <p>Read the quickstart tutorial to learn how to get set up with the required software and tools, and how to launch a nf-core pipeline.</p>
      <a class="btn btn-lg btn-success arrow-hover" href="/usage/introduction"><span><i class="fad fa-graduation-cap me-1"></i> Quickstart Tutorial</span></a>
    </div>
    <div class="col-md-4 px-md-5 mb-5 mb-md-0">
      <h3>Get into the code</h3>
      <p>If you're interested in contributing to sanger-tol, take a look at the developer documentation to see what's required.</p>
      <a class="btn btn-lg btn-success arrow-hover" href="/docs/contributing/adding_pipelines"><span><i class="fad fa-code me-1"></i> Developer docs</span></a>
    </div>
  </div>
</div>

<!--
<div id="community" class="homepage-usedby">
  <div class="container py-5">
    <h2>
      <a class="btn btn-success float-end d-none d-md-inline" href="/community#organisations">See a complete list &raquo;</a>
      <a href="/community#organisations">Used by groups all over the world</a>
    </h2>
    <p>The nf-core community is spread all over the globe and includes a large
      number of contributing users.</p>
    <p><a class="btn btn-success d-inline d-md-none" href="/community#organisations">See a complete list &raquo;</a></p>
    <div class="homepage_contrib_logos">
      <?php foreach (array_slice($contributors_img_list, 0, 8) as $img) {
          echo $img;
      } ?>
    </div>
  </div>
</div>
-->

<script type="text/javascript">
// List of remaining contributor logos which have not yet been shown
var contributors_imgs = <?php echo json_encode(array_slice($contributors_img_list, 8)); ?>;

<?php // Javascript for moment time zone support

if ($event['start_time']) {
    echo '
    $("[data-timestamp]").each(function(){
      var timestamp = $(this).data("timestamp");
      var local_time = moment.tz(timestamp, "X", moment.tz.guess());
      $(this).text(local_time.format("HH:mm z, LL"));
    });
    ';
} ?>
</script>

<?php
$md_github_url = false;
include '../includes/footer.php';


?>
