<?php

// Check if we have $_GET['t'] - ie, URL was /join/something
$config = parse_ini_file("../config.ini");
$join_redirects = [
  'slack' => $config['slack_invite_url']
];
if(isset($_GET['t'])){
  $redirect = rtrim(str_replace('join/', '', $_GET['t']), '/');
  if(array_key_exists($redirect, $join_redirects)){
    header('Location: '.$join_redirects[$redirect]);
  } else {
    header('HTTP/1.1 404 Not Found');
    include('404.php');
    die();
  }
}

$title = 'Join nf-core';
$subtitle = 'Read about the different ways you can get involved with nf-core';
include('../includes/header.php');
?>

<p>We use a few different tools to organise the nf-core community -
    you are welcome to join us at any or all!</p>
<p class="text-info small">
  <i class="fas fa-exclamation-triangle mr-1"></i>
  All nf-core community members are expected to adhere to the nf-core <a href="/code_of_conduct" class="text-info link-underline">code of conduct</a>
</p>
<p class="small" style="color: #D93F66;">
  <i class="fab fa-gitter mr-1"></i>
  If your question is about Nextflow and not directly related to nf-core, please use the main <a style="color: #D93F66;" class="link-underline" href="https://gitter.im/nextflow-io/nextflow">Nextflow Gitter chat</a>.
</p>
<p class="text-center">
  <a href="https://nf-co.re/join/slack" class="mb-2 btn btn-success">
    <i class="fab fa-slack"></i> Join Slack
  </a>
  <a href="https://github.com/nf-core" class="mb-2 btn btn-dark">
    <i class="fab fa-github"></i> nf-core on GitHub
  </a>
  <a href="https://twitter.com/nf_core" class="mb-2 btn btn-info" style="background-color: #4AA1EB;">
    <i class="fab fa-twitter"></i> @nf_core on twitter
  </a>
  <a href="https://www.youtube.com/c/nfcore" class="mb-2 btn btn-danger" style="background-color: #EA3C1E;">
    <i class="fab fa-youtube"></i> YouTube
  </a>
  <a href="https://groups.google.com/forum/#!forum/nf-core" class="mb-2 btn btn-info" style="background-color: #5CA5EF;">
    <i class="fas fa-envelope"></i> Google Groups email list
  </a>
</p>

<h1>
  <a class="text-success text-decoration-none" href="https://nfcore.slack.com/" target="_blank">
    <span class="join-pull-img"><img height="120px" src="/assets/img/slack.svg" /></span>
    Slack
  </a>
</h1>
<p class="text-info small">
  <i class="fas fa-question-circle mr-1"></i>
  If you would like help with running nf-core pipelines, Slack is the best place to start.
</p>
<p>Slack is a real-time messaging tool, with discussion split into channels and groups.
We use it to provide help to people running nf-core pipelines, as well as discussing development ideas.
You can join the nf-core slack <a href="https://nf-co.re/join/slack">here</a>.</p>
<p>Once you have registered, you can access the nf-core slack at <a href="https://nfcore.slack.com/">https://nfcore.slack.com/</a>
  <em class="small text-muted">(NB: No hyphen!)</em></p>

<h1>
  <a class="text-success text-decoration-none" href="https://github.com/nf-core/" target="_blank">
    <span class="join-pull-img join-gh-logo"><img height="120px" src="/assets/img/github.svg" /></span>
    GitHub organisation
  </a>
</h1>
<p class="text-info small">
  <i class="fas fa-bug mr-1"></i>
  If you encounter a bug or have a suggestion, please create an issue on the repository for that pipeline.
</p>
<p>We use GitHub to manage all of the code written for nf-core.
It's a fantastic platform and provides a huge number of tools.
We have a GitHub organisation called <a href="https://github.com/nf-core/">nf-core</a> which anyone can join:
drop us a note <a href="https://github.com/nf-core/nf-co.re/issues/3">here</a> or anywhere and we'll send you an invite.
</p>

<h1>
  <a class="text-success text-decoration-none" href="https://twitter.com/nf_core" target="_blank">
    <span class="join-pull-img"><img height="120px" src="/assets/img/twitter.svg" /></span>
    Twitter
  </a>
</h1>
<p>The <a href="https://twitter.com/nf_core">@nf_core</a> is a low-volume account that sends
    automated tweets whenever a new pipeline release is tagged. Relevant news and events are occasionally also tweeted.
    See <a href="https://twitter.com/nf_core">https://twitter.com/nf_core</a></p>

<h1>
  <a class="text-success text-decoration-none" href="https://www.youtube.com/c/nfcore" target="_blank">
    <span class="join-pull-img"><img width="144px" src="/assets/img/youtube.svg" /></span>
    YouTube
  </a>
</h1>
<p>The <a href="https://www.youtube.com/c/nfcore">nf-core</a> YouTube account is used for tutorial videos
    and has a playlist to collect recordings of presentations about nf-core from across the web.
    See <a href="https://www.youtube.com/c/nfcore">https://www.youtube.com/c/nfcore</a>
</p>

<h1>
  <a class="text-success text-decoration-none" href="https://groups.google.com/forum/#!forum/nf-core" target="_blank">
    <img height="120px" src="/assets/img/google_groups.svg" class="join-pull-img" />
    Google Groups
  </a>
</h1>
<p>If Slack isn't your thing and you prefer traditional email lists, head over to
    <a href="https://groups.google.com/forum/#!forum/nf-core">https://groups.google.com/forum/#!forum/nf-core</a> and send us a message!</p>

<?php include('../includes/footer.php');
