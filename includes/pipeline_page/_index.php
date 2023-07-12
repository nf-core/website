<?php

require_once '../includes/functions.php';
usort($pipeline->releases, 'rsort_releases');

########
## Redirect old style URLs
########
if (count($path_parts) > 2 && $path_parts[1] == 'docs' && $path_parts[2] == 'usage') {
    header('Location: /' . $pipeline->name . '/usage');
}
if (count($path_parts) > 2 && $path_parts[1] == 'docs' && $path_parts[2] == 'output') {
    header('Location: /' . $pipeline->name . '/output');
}

########
## Configure page header
########
$title = 'sanger-tol/<br class="d-sm-none">' . $pipeline->name;
$subtitle = $pipeline->description;
$content = '';
$schema_content = '';
$import_chartjs = true;
$no_auto_toc = true;

$subpage_names = [
    '' => 'Introduction',
    'results' => 'Results',
    'usage' => 'Usage docs',
    'parameters' => 'Parameters',
    'output' => 'Output docs',
    'releases_stats' => 'Releases & Statistics',
];

# Header - keywords
$header_html = '<p class="mb-5">';
foreach ($pipeline->topics as $keyword) {
    $header_html .= '<a href="/pipelines?q=' . $keyword . '" class="badge pipeline-topic">' . $keyword . '</a> ';
}
$header_html .= '</p>';

// Highlight any search terms if we have them
if (isset($_GET['q']) && strlen($_GET['q'])) {
    $title = preg_replace('/(' . $_GET['q'] . ')/i', "<mark>$1</mark>", $title);
    $subtitle = preg_replace('/(' . $_GET['q'] . ')/i', "<mark>$1</mark>", $subtitle);
    $header_html = preg_replace('/(' . $_GET['q'] . ')/i', "<mark>$1</mark>", $header_html);
}

########
## Work out the pipeline release
########

# Set defaults (Readme tab)
$pagetab = ''; # empty string is home / readme
$release = 'dev';
$release_hash = null;
$latest_release = 'dev';
if (count($pipeline->releases) > 0) {
    $release = $latest_release = $pipeline->releases[0]->tag_name;
    $release_url = $pipeline->releases[0]->html_url;
    $release_hash = $pipeline->releases[0]->tag_sha;
}
# Find release from URL if set
if (count($path_parts) > 1) {
    foreach ($pipeline->releases as $r) {
        if ($path_parts[1] === $r->tag_name) {
            $release = $r->tag_name;
            $release_url = $releases->html_url;
            $release_hash = $r->tag_sha;
        }
    }
    if ($path_parts[1] == 'dev') {
        $release = 'dev';
        $release_url = null;
        $release_hash = null;
    }
}
########
## Load stats from database
########

$config = parse_ini_file('../config.ini');
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// get stats for current pipeline
$sql = "SELECT * FROM nfcore_pipelines WHERE name = '" . $pipeline->name . "'";
$pipeline_metrics = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipeline_metrics = mysqli_fetch_all($result, MYSQLI_ASSOC)[0];
        // Free result set
        mysqli_free_result($result);
    }
}
// get traffic stats for current pipeline
$sql =
    "SELECT * FROM github_traffic_stats WHERE pipeline_id = '" . $pipeline_metrics['id'] . "' ORDER BY timestamp DESC";
$traffic_stats = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $traffic_stats = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}

// get contributor stats for current pipeline
$sql =
    "SELECT * FROM github_pipeline_contrib_stats WHERE pipeline_id = '" .
    $pipeline_metrics['id'] .
    "' ORDER BY week_date DESC";
$contributor_stats = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $contributor_stats = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    }
}

mysqli_close($conn);
########
## Load and cache the pipeline JSON schema if we have one
########
# Make cache file names
$gh_pipeline_schema_fn = dirname(dirname(__FILE__)) . "/api_cache/json_schema/{$pipeline->name}/{$release}.json";
$gh_pipeline_no_schema_fn =
    dirname(dirname(__FILE__)) . "/api_cache/json_schema/{$pipeline->name}/{$release}.NO_SCHEMA";
# Build directories if needed
if (!is_dir(dirname($gh_pipeline_schema_fn))) {
    mkdir(dirname($gh_pipeline_schema_fn), 0777, true);
}
// Try to fetch the nextflow_schema.json file for the selected release, if not already cached or release==dev. Decides later if Launch button is included on the page or not.
if ((!file_exists($gh_pipeline_schema_fn) && !file_exists($gh_pipeline_no_schema_fn)) || $release == 'dev') {
    $api_opts = stream_context_create(['http' => ['method' => 'GET', 'header' => ['User-Agent: PHP']]]);
    $gh_launch_schema_url = "https://api.github.com/repos/sanger-tol/{$pipeline->name}/contents/nextflow_schema.json?ref={$release}";
    $gh_launch_schema_json = file_get_contents($gh_launch_schema_url, false, $api_opts);
    if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
        echo '<script>console.log("Sent request to ' .
            $gh_launch_schema_url .
            '"," got http response header:",' .
            json_encode($http_response_header, JSON_HEX_TAG) .
            ')</script>';
        # Remember for next time
        if (strpos($http_response_header[0], 'HTTP/1.1 404') !== false) {
            file_put_contents($gh_pipeline_no_schema_fn, '');
        }
    } else {
        # Parse out file content
        $gh_launch_schema_response = json_decode($gh_launch_schema_json, true);
        $gh_launch_schema = base64_decode($gh_launch_schema_response['content'], true);
        # Save cache
        file_put_contents($gh_pipeline_schema_fn, $gh_launch_schema);
    }
}

########
## Figure out what page we're rendering
########

# URL path to readme - redirect to pipeline root
if (endswith($_GET['path'], '/README')) {
    header('Location: /' . substr($_GET['path'], 0, -7));
    exit();
}

# Usage docs
if (endswith($_GET['path'], '/usage')) {
    $pagetab = 'usage';
    $filename = 'docs/usage.md';
    $md_trim_before = '# Introduction';
}
# Usage docs
elseif (endswith($_GET['path'], '/parameters')) {
    $pagetab = 'parameters';
    require_once 'docs_schema.php';
}
# Output docs
elseif (endswith($_GET['path'], '/output')) {
    $pagetab = 'output';
    $filename = 'docs/output.md';
    $md_trim_before = '# Introduction';
}
# Example output
elseif (endswith($_GET['path'], '/results')) {
    $pagetab = 'results';
    require_once 'results.php';
}
# Releases + Stats
elseif (endswith($_GET['path'], '/releases_stats')) {
    $pagetab = 'releases_stats';
    require_once 'pipeline_releases_stats.php';
}
# Some other URL pattern that we don't recognise - 404
elseif ($_GET['path'] != $pipeline->name && $_GET['path'] != $pipeline->name . '/' . $release) {
    $protocol = isset($_SERVER['HTTPS']) && !empty($_SERVER['HTTPS']) ? 'https://' : 'http://';
    $url_string = trim(str_replace($pipeline->name, '', $_GET['path']), '/');
    $url_string = trim(str_replace($release, '', $url_string), '/');
    header('HTTP/1.1 404 Not Found');
    $suggestion_404_urls = [
        $protocol . $_SERVER['HTTP_HOST'] . '/' . $pipeline->name,
        'https://github.com/sanger-tol/' . $pipeline->name . '/blob/' . $release . '/' . $url_string,
        'https://github.com/snger-tol/' . $pipeline->name . '/blob/' . $release . '/' . $url_string . '.md',
    ];
    include '404.php';
    die();
}
# Homepage,
if ($pagetab == '') {
    $filename = 'README.md';
    $md_trim_before = '# Introduction';
}

# Prep local cache and variables for docs markdown
// NB: $content rendered in header.php
if (in_array($pagetab, ['', 'usage', 'output'])) {
    require_once 'docs_md.php';
}

# Main page nav and header
$url_base = '/' . $pipeline->name . '/' . $release;
$no_print_content = true;
$mainpage_container = false;
include '../includes/header.php';

########
# Get details for the Call To Action button
########
# Text, URL and icon for the CTA button
$cta_txt = $release == 'dev' ? 'See the latest code' : 'See version ' . $release;
$cta_url = $pipeline->html_url;
$cta_icon = $release == 'dev' ? '<i class="fad fa-construction me-1"></i> ' : '<i class="fas fa-tags me-1"></i> ';
if (file_exists($gh_pipeline_schema_fn)) {
    $cta_txt = $release == 'dev' ? 'Launch development version' : 'Launch version ' . $release;
    $cta_url = '/launch?pipeline=' . $pipeline->name . '&release=' . $release;
    $cta_icon = '<i class="fad fa-rocket-launch me-1"></i> ';
}
# Build button
$cta_btn = '<a href="' . $cta_url . '" class="btn btn-success btn-lg d-print-none">' . $cta_icon . $cta_txt . '</a>';

########
# Warning alert box
########
$pipeline_warning = '';
$latest_stable =
    'The latest stable release is <a href="/' .
    $pipeline->name .
    '/' .
    $latest_release .
    '"><code>v' .
    $latest_release .
    '</code></a>';
if (count($pipeline->releases) == 0) {
    $pipeline_warning =
        '<div class="alert alert-danger">This pipeline is currently in development and does not yet have any stable releases.</div>';
} elseif ($release == 'dev') {
    $pipeline_warning =
        '<div class="alert alert-warning">You are viewing the development version pages for this pipeline. ' .
        $latest_stable .
        '</div>';
} elseif ($release != $latest_release) {
    $pipeline_warning =
        '<div class="alert alert-warning">These pages are for an old version of the pipeline (<code>v' .
        $release .
        '</code>). ' .
        $latest_stable .
        '</div>';
}
if ($pipeline->archived) {
    $pipeline_warning =
        '<div class="alert alert-warning">This pipeline has been archived and is no longer being actively maintained.</div>';
}

########
# Extra HTML for the header - tags and GitHub URL
########
?>

<div class="container-fluid mainpage-subheader-heading chevron-down">
  <div class="container text-center">
    <?php echo $pipeline_warning; ?>
    <p><?php echo $cta_btn; ?></p>
    <p><a href="<?php echo $pipeline->html_url; ?>" class="subheader-link"><i class="fab fa-github"></i> <?php echo $pipeline->html_url; ?></a></p>
  </div>
</div>

<div class="container-xxl main-content">

    <ul class="nav nav-fill nfcore-subnav justify-content-start justify-content-md-around d-print-none">
    <li class="nav-item dropdown d-block d-lg-none" style="z-index:10000;">
        <button class="btn btn-secondary dropdown-toggle" type="button" data-bs-toggle="dropdown" aria-expanded="false">
            <?php echo $subpage_names[$pagetab]; ?>
        </button>
        <div class="dropdown-menu" aria-labelledby="dropdownMenuButton">
        <a class="dropdown-item<?php if ($pagetab == '') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>"><i class="fas fa-sign-in me-1"></i>
        <?php echo $subpage_names['']; ?></a>
        <a class="dropdown-item<?php if ($pagetab == 'results') {
            echo ' active';
        } ?>" href="/<?php echo $pipeline->name; ?>/results"><i class="f-regular fa-folder-tree fa-lg me-1"></i>
            <?php echo $subpage_names['results']; ?>
        </a>
        <a class="dropdown-item<?php if ($pagetab == 'usage') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>/usage"><i class="far fa-book me-1"></i>
            <?php echo $subpage_names['usage']; ?>
        </a>
        <?php if (file_exists($gh_pipeline_schema_fn)): ?>
        <a class="dropdown-item<?php if ($pagetab == 'parameters') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>/parameters"><i class="far fa-book me-1"></i>
            <?php echo $subpage_names['parameters']; ?>
        </a>
        <?php endif; ?>
        <a class="dropdown-item<?php if ($pagetab == 'output') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>/output"><i class="far fa-book me-1"></i>
            <?php echo $subpage_names['output']; ?>
        </a>
        <a class="dropdown-item<?php if ($pagetab == 'releases_stats') {
            echo ' active';
        } ?>" href="/<?php echo $pipeline->name; ?>/releases_stats"><i class="fas fa-chart-line me-1"></i>
        <?php echo $subpage_names['releases_stats']; ?><span class="d-none d-sm-inline">istic</span>s</a>
    </li>


    <li class="nav-item d-none d-lg-block">
        <a class="nav-link<?php if ($pagetab == '') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>"><i class="fas fa-sign-in me-1"></i> <?php echo $subpage_names['']; ?></a>
    </li>
    <?php if (isset($release_hash) && $release_hash): ?>
        <li class="nav-item d-none d-lg-block">
        <a class="nav-link<?php if ($pagetab == 'results') {
            echo ' active';
        } ?>" href="/<?php echo $pipeline->name; ?>/results"><i class="fa-regular fa-folder-tree fa-lg me-1"></i>
        <?php echo $subpage_names['results']; ?></span></a>
        </li>
    <?php endif; ?>
    <li class="nav-item d-none d-lg-block">
      <a class="nav-link<?php if ($pagetab == 'usage') {
          echo ' active';
      } ?>" href="<?php echo $url_base; ?>/usage"><i class="far fa-book me-1"></i>
      <?php echo $subpage_names['usage']; ?></a>
    </li>
    <?php if (file_exists($gh_pipeline_schema_fn)): ?>
      <li class="nav-item d-none d-lg-block">
        <a class="nav-link<?php if ($pagetab == 'parameters') {
            echo ' active';
        } ?>" href="<?php echo $url_base; ?>/parameters"><i class="far fa-book me-1"></i>
        <?php echo $subpage_names['parameters']; ?></a>
      </li>
    <?php endif; ?>
    <li class="nav-item d-none d-lg-block">
      <a class="nav-link<?php if ($pagetab == 'output') {
          echo ' active';
      } ?>" href="<?php echo $url_base; ?>/output"><i class="far fa-book me-1"></i>
      <?php echo $subpage_names['output']; ?></a>
    </li>

    <li class="nav-item d-none d-lg-block">
      <a class="nav-link<?php if ($pagetab == 'releases_stats') {
          echo ' active';
      } ?>" href="/<?php echo $pipeline->name; ?>/releases_stats"><i class="fas fa-chart-line me-1"></i>
      <?php echo $subpage_names['releases_stats']; ?></a>
    </li>

    <?php if ($pagetab != 'releases_stats'): ?>
      <li class="nav-item pt-1 ps-3">
        <div class="input-group input-group-sm">
          <label class="input-group-text" for="version_select"><i class="fas fa-tags"></i></label>
          <select class="form-select form-select-sm" id="version_select" data-pipeline="<?php echo $pipeline->name; ?>">
            <?php
            $releases = [];
            foreach ($pipeline->releases as $r) {
                $releases[$r->tag_name] = $r->tag_sha;
            }
            $releases['dev'] = '';
            foreach ($releases as $r => $h) {
                $selected = $r === $release ? 'selected="selected"' : '';
                echo '<option value="/' .
                    $pipeline->name .
                    '/' .
                    $r .
                    '/' .
                    $pagetab .
                    '" ' .
                    $selected .
                    '>' .
                    $r .
                    '</option>';
            }
            ?>
          </select>
        </div>
      </li>
    <?php endif; ?>

  </ul>

  <?php
  ########
  # Make a row with a column for content
  ########
  if (in_array($pagetab, ['results'])) {
      echo '<div class="row ms-lg-5"><div class="col-12">';
  } elseif (in_array($pagetab, ['', 'releases_stats'])) {
      echo '<div class="row ms-lg-5"><div class="col-12 col-lg-9">';
  } else {
      echo '<div class="row ms-lg-5 flex-wrap-reverse flex-lg-wrap"><div class="col-12 col-lg-9">';
  }

  ########
  # Print content
  ########
  # Add on the rendered schema docs (empty string if we don't have it)
  if (preg_match('/<!-- params-docs -->/', $content)) {
      $content = preg_replace('/<!-- params-docs -->/', $schema_content, $content);
  } else {
      $content .= $schema_content;
  }

  echo '<div class="rendered-markdown pipeline-page-content">' . $content . '</div>';

  echo '</div>'; # end of the content div
  if (in_array($pagetab, ['results'])) {
      echo '<div><div">';
  } elseif (in_array($pagetab, ['', 'releases_stats'])) {
      echo '<div class="col-12 col-lg-3 ps-2 h-100 sticky-md-top"><div class="side-sub-subnav">';
  } else {
      echo '<div class="col-12 col-lg-3 ps-2 h-100 sticky-top"><div class="side-sub-subnav">';
  }

  # Pipeline homepage & releases - key stats
  if (in_array($pagetab, ['', 'releases_stats'])) {
      require_once 'sidebar.php';
  }
  # Documentation - ToC
  elseif (in_array($pagetab, ['usage', 'output', 'parameters'])) {
      $toc .= '<nav class="toc mt-2 auto-toc border-start d-print-none">';
      $toc .= generate_toc($content);
      $toc_md = '';
      # Add on the action buttons for the parameters docs
      if ($pagetab == 'parameters') {
          $toc .= '
            <div class="btn-group w-100 mt-2 mb-1 ms-1 d-none d-lg-block" role="group">
            <button class="btn btn-sm btn-outline-secondary" data-bs-toggle="collapse" data-bs-target=".param-docs-help-text" aria-expanded="false">
                <i class="fas fa-question-circle me-1"></i> Show all help
            </button>
            <button class="btn btn-sm btn-outline-secondary btn-show-hidden-params" data-bs-toggle="collapse" data-bs-target=".param-docs-hidden, .toc .collapse, .hidden_params_alert" aria-expanded="false">
                <span class="collapse show"><i class="fas fa-eye-slash"></i> Show hidden params</span>
                <span class="collapse"><i class="fas fa-eye"></i> Hide hidden params</span>
            </button>
            </div>';
          $toc_md = '
            <a class="dropdown-item" data-bs-toggle="collapse" data-bs-target=".param-docs-help-text" aria-expanded="false">
                <i class="fas fa-question-circle me-1"></i> Show all help
            </a>
            <a class="dropdown-item btn-show-hidden-params" data-bs-toggle="collapse" data-bs-target=".param-docs-hidden, .toc .collapse, .hidden_params_alert" aria-expanded="false">
                <span class="collapse show"><i class="fas fa-eye-slash"></i> Show hidden params</span>
                <span class="collapse"><i class="fas fa-eye"></i> Hide hidden params</span>
            </a>';
      }
      if ($pagetab == 'output') {
          $toc .= '
            <div class="btn-group w-100 mt-2 mb-1 ms-1 d-none d-lg-block" role="group">
            <div class="text-center ms-1">
                <button class="btn btn-sm btn-outline-secondary expand-details w-100">
                <span class=""><i class="fas fa-arrows-v"></i> Expand all output file details</span>
                </button>
            </div>
            </div>';
          $toc_md .= '
                <a class=" dropdown-item expand-details">
                <span class=""><i class="fas fa-arrows-v"></i> Expand all output file details</span>
                </a>';
      }
      $toc = preg_replace('/<!-- tock_md_button_placeholder -->/', $toc_md, $toc);
      # Back to top link
      $toc .=
          '<p class="small text-end mt-3 d-none d-lg-block"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
      $toc .= '</nav>';
      echo $toc;
  }
  echo '</div></div>'; # end of the sidebar col
  echo '</div>'; # end of the row
  include '../includes/footer.php';

