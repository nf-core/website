<?php

// Check if we have $_GET['t'] - ie, URL was /join/something
$config = parse_ini_file("../config.ini");
$join_redirects = [
  'slack' => $config['slack_invite_url']
];
if(isset($_GET['t'])){
  $redirect = str_replace('join/', '', $_GET['t']);
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
  <i class="fas fa-question-circle"></i>
  If you would like help with running nf-core pipelines, Slack is the best place to start.
</p>
<p class="text-info small">
  <i class="fas fa-bug"></i>
  If you encounter a bug or have a suggestion, please create an issue on the repository for that pipeline.
</p>
<p class="text-center">
  <a href="https://nf-co.re/join/slack" class="mb-2 btn btn-success">
    <i class="fab fa-slack"></i> Get Slack invite
  </a>
  <a href="https://groups.google.com/forum/#!forum/nf-core" class="mb-2 btn btn-info" style="background-color: #5CA5EF;">
    <i class="fas fa-envelope"></i> Google Groups email list
  </a>
  <a href="https://github.com/nf-core" class="mb-2 btn btn-dark">
    <i class="fab fa-github"></i> nf-core on GitHub
  </a>
  <a href="https://twitter.com/nf_core" class="mb-2 btn btn-info" style="background-color: #4AA1EB;">
    <i class="fab fa-twitter"></i> @nf_core on twitter
  </a>
</p>

<h1>
  <img height="120px" src="/assets/img/slack.svg" class="float-right bg-white pl-4" />
  Slack
</h1>
<p>Slack is a real-time messaging tool, with discussion split into channels and groups.
We use it to provide help to people running nf-core pipelines, as well as discussing development ideas.
You can join the nf-core slack by getting an invite <a href="https://nf-co.re/join/slack">here</a>.</p>
<p>Once you have registered, you can access the nf-core slack at <a href="https://nf-core.slack.com/">https://nf-core.slack.com/</a></p>

<h1>
  <img height="120px" src="/assets/img/github.svg" class="float-right bg-white pl-4" />
  GitHub organisation
</h1>
<p>We use GitHub to manage all of the code written for nf-core.
It's a fantastic platform and provides a huge number of tools.
We have a GitHub organisation called <a href="https://github.com/nf-core/">nf-core</a> which anyone can join:
drop us a note <a href="https://github.com/nf-core/nf-co.re/issues/3">here</a> or anywhere and we'll send you an invite.
</p>

<h1>
  <img height="120px" src="/assets/img/twitter.svg" class="float-right bg-white pl-4" />
  Twitter
</h1>
<p>The <a href="https://twitter.com/nf_core">@nf_core</a> is a low-volume account that sends
    automated tweets whenever a new pipeline release is tagged. Relevant news and events are occasionally also tweeted.
    See <a href="https://twitter.com/nf_core">https://twitter.com/nf_core</a></p>

<h1>
  <img height="120px" src="/assets/img/google_groups.svg" class="float-right bg-white pl-4" />
  Google Groups
</h1>
<p>If slack isn't your thing and you prefer traditional email lists, head over to
    <a href="https://groups.google.com/forum/#!forum/nf-core">https://groups.google.com/forum/#!forum/nf-core</a> and send us a message!</p>

<?php include('../includes/footer.php');
