<?php
$title = 'Events';
$subtitle = 'Details of past and future nf-core meetups.';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/nf-core-events.yaml';
include('../includes/header.php');

require_once("../Spyc.php");
$events = spyc_load_file('../nf-core-events.yaml');

include('../includes/footer.php');
