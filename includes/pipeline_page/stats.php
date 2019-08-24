<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.

ob_start();

echo '<p class="lead small mt-5">Coming soon..</p>';

$content = ob_get_contents();
ob_end_clean();
