<?php
$title = 'Code of conduct';
$subtitle = 'The nf-core code of conduct';
$markdown_fn = '../CODE_OF_CONDUCT.md';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/CODE_OF_CONDUCT.md';
$no_print_content = true;
include('../includes/header.php');
?>

The nf-core project follows the <a href="https://www.contributor-covenant.org/" target="_blank">Contributor Covenant</a>,
which can be seen below.

<?php
echo $content;
include('../includes/footer.php');
