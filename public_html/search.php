<?php

$title = 'Search nf-core';
$search_term = $_GET['q'];
$search_term = htmlspecialchars($search_term, ENT_QUOTES, 'UTF-8');
$subtitle = 'Searching for the term <code>' . $search_term . '</code>';
$mainpage_container = false;
include '../includes/header.php';

?>
<div class="mainpage-subheader-heading chevron-down">
  <div class="container  my-3">
    <form class="searchbar_form" action="search" method="get">
      <div class="col-7 col-lg-3 mx-auto">
        <div class="input-group">
          <input type="search" class="form-control" placeholder="Search" name="q" required value="<?php echo $search_term; ?>">
          <button type="submit" class="btn btn-outline-success">Search</button>
        </div>
      </div>
    </form>
  </div>
</div>


<div class="container main-content">

  <?php
  include '../includes/search_results.php';

  # DEBUG
  # echo '<pre>'.print_r($search_results, true).'</pre>';

  if (count($search_results['pipelines']) == 0 && count($search_results['documentation']) == 0) {
      echo '<h1>No results</h1>';
      echo '<p>Sorry, no search results were found. Please try another term.</p>';
  }

  if (count($search_results['pipelines']) > 0) {
      echo '<h1>Pipelines</h1>';
      foreach ($search_results['pipelines'] as $result) {
          $wf = $result['pipeline']; ?>
      <div class="card search-page-result mb-2">
        <div class="card-body">
          <h5 class="card-title mb-0"><a href="<?php echo $wf->url .
              "?q=$search_term"; ?>"><?php echo $wf->full_name; ?></a></h5>
          <?php if (count($wf->topics) > 0): ?>
            <p class="topics mb-0">
              <?php foreach ($wf->topics as $topic): ?>
                <a href="/pipelines?q=<?php echo $topic; ?>" class="badge pipeline-topic"><?php echo $topic; ?></a>
              <?php endforeach; ?>
            </p>
          <?php endif; ?>
          <p class="card-text text-muted small"><?php echo $wf->description; ?></p>
        </div>
      </div>
    <?php
      }
  }

  if (count($search_results['documentation']) > 0) {
      echo '<h1>Documentation</h1>';
      foreach ($search_results['documentation'] as $result) { ?>
      <div class="card search-page-result mb-2">
        <div class="card-body">
          <h5 class="card-title mb-0"><a href="<?php echo $result['url'] . "?q=$search_term"; ?>"><?php echo $result[
    'title'
]; ?></a></h5>
          <h6 class="text-muted"><?php echo $result['subtitle']; ?></h6>
          <p class="card-text text-muted small"><?php echo $result['match_string']; ?></p>
        </div>
      </div>
  <?php }
  }

  include '../includes/footer.php';

