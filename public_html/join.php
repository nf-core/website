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
<p class="text-center">
  <a href="https://nf-co.re/join/slack" class="mb-2 btn btn-lg btn-success">
    <i class="fab fa-slack"></i> Join nf-core Slack
  </a>
  <a href="https://github.com/nf-core" class="mb-2 btn btn-lg btn-dark">
    <i class="fab fa-github"></i> nf-core on GitHub
  </a>
  <a href="https://twitter.com/nf_core" class="mb-2 btn btn-lg btn-info" style="background-color: #4AA1EB;">
    <i class="fab fa-twitter"></i> @nf_core on twitter
  </a>
  <a href="https://www.youtube.com/c/nf-core" class="mb-2 btn btn-lg btn-danger" style="background-color: #EA3C1E;">
    <i class="fab fa-youtube"></i> nf-core on YouTube
  </a>
</p>
<p class="text-info small">
  <i class="fas fa-exclamation-triangle mr-1"></i>
  All nf-core community members are expected to adhere to the nf-core <a href="/code_of_conduct" class="text-info link-underline">code of conduct</a>
</p>
<div class="alert alert-info text-center">
  <i class="fab fa-gitter mr-1"></i>
  If your question is about Nextflow and not directly related to nf-core, please use the main <a class="link-underline" href="https://gitter.im/nextflow-io/nextflow">Nextflow Gitter chat</a>.
</div>

<h1 id="slack" class="mt-5"><a href="#slack" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
  <span class="join-pull-img"><img height="120px" src="/assets/img/slack.svg" /></span>
  Slack
</h1>
<p class="text-info small">
  <i class="fas fa-question-circle mr-1"></i>
  If you would like help with running nf-core pipelines, Slack is the best place to start.
</p>
<p>Slack is a real-time messaging tool, with discussion split into channels and groups.
We use it to provide help to people running nf-core pipelines, as well as discussing development ideas.</p>
<p>Once you have registered, you can access the nf-core slack at <a href="https://nfcore.slack.com/">https://nfcore.slack.com/</a>
  <em class="small text-muted">(NB: No hyphen!)</em></p>
<p><a class="btn btn-success mr-2" href="https://nf-co.re/join/slack"><i class="fab fa-slack"></i> Get invite to nf-core Slack</a>
</p>

<h1 id="github" class="mt-5"><a href="#github" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
  <span class="join-pull-img join-gh-logo"><img height="120px" src="/assets/img/github.svg" /></span>
  GitHub organisation
</h1>
<p class="text-info small">
  <i class="fas fa-bug mr-1"></i>
  If you encounter a bug or have a suggestion, please create an issue on the repository for that pipeline.
</p>
<p>We use GitHub to manage all of the code written for nf-core.
It's a fantastic platform and provides a huge number of tools.
We have a GitHub organisation called <a href="https://github.com/nf-core/">nf-core</a> which we use to give access to
collaborators for editing and reviewing code.</p>
<p>To join the GitHub organisation, please drop us a note <a href="https://github.com/nf-core/nf-co.re/issues/3">here</a>
or on the <a href="https://nfcore.slack.com/channels/github-invitations"><code>#github-invitations</code> Slack channel</a> and we'll send an invite.</p>

<p><a href="https://github.com/nf-core" class="mb-2 btn btn-dark">
  <i class="fab fa-github"></i> nf-core on GitHub
</a></p>

<h1 id="twitter" class="mt-5"><a href="#twitter" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
  <span class="join-pull-img"><img height="120px" src="/assets/img/twitter.svg" /></span>
  Twitter
</h1>
<p>The <a href="https://twitter.com/nf_core">@nf_core</a> twitter account sends automated tweets whenever
  a new pipeline release is tagged. Relevant news and events are also tweeted.</p>
<a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-size="large">Follow @nf_core</a>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

<h1 id="youtube" class="mt-5"><a href="#youtube" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
  <span class="join-pull-img"><img width="120px" src="/assets/img/youtube.svg" /></span>
  YouTube
</h1>
<p>The <a href="https://www.youtube.com/c/nf-core">nf-core</a> YouTube channel is used for tutorial videos
    and has a playlist to collect recordings of presentations about nf-core from across the web.
    See <a href="https://www.youtube.com/c/nf-core">https://www.youtube.com/c/nf-core</a>
</p>
<script src="https://apis.google.com/js/platform.js"></script>
<div class="g-ytsubscribe" data-channelid="UCfdY8g6eA6IuYZ8oyyNsabQ" data-layout="full" data-count="default"></div>

<?php include('../includes/footer.php');
