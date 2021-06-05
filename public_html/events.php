<?php

# To get parse_md_front_matter() and sanitise_date_meta() functions
require_once('../includes/functions.php');

$md_base = dirname(dirname(__file__)) . "/markdown/";
$event_type_classes = array(
  'hackathon' => 'primary',
  'talk' => 'success',
  'poster' => 'secondary',
  'tutorial' => 'info',
  'workshop' => 'light'
);
$event_type_icons = array(
  'hackathon' => 'fas fa-laptop-code',
  'talk' => 'fas fa-presentation',
  'poster' => 'far fa-image',
  'tutorial' => 'fas fa-graduation-cap',
  'workshop' => 'fas fa-chalkboard-teacher'
);

function print_events($events, $is_past_event)
{
  global $event_type_classes;
  global $event_type_icons;

  foreach ($events as $idx => $event) :
    # Nice date strings
    $date_string = date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) . ' - ' . date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    if (date('mY', $event['start_ts']) == date('mY', $event['end_ts'])) {
      $date_string = date('j<\s\u\p>S</\s\u\p> ', $event['start_ts']) . ' - ' . date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    }
    if (date('dmY', $event['start_ts']) == date('dmY', $event['end_ts'])) {
      $date_string = date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    }
    $colour_class = $event_type_classes[strtolower($event['type'])];
    $icon_class = $event_type_icons[strtolower($event['type'])];
?>

    <!-- Event Card -->
    <div class="card my-4 border-top-0 border-right-0 border-bottom-0 border-<?php echo $colour_class ?>">
      <div class="card-body <?php if ($is_past_event) {
                              echo 'py-2';
                            } ?>">
        <h5 class="my-0 py-0">
          <a class="text-success" href="<?php echo $event['url']; ?>"><?php echo $event['title']; ?></a>
          <small><span class="badge badge-<?php echo $colour_class ?> float-right small"><i class="<?php echo $icon_class ?> mr-1"></i><?php echo ucfirst($event['type']); ?></span></small>
        </h5>
        <?php if (array_key_exists('subtitle', $event)) {
          $tm = $is_past_event ? 'text-muted' : '';
          echo '<p class="mb-0 ' . $tm . '">' . $event['subtitle'] . '</p>';
        }
        if (!$is_past_event) : ?>
          <h6 class="small text-muted"><?php echo $date_string; ?></h6>

          <?php if (array_key_exists('description', $event)) {
            echo '<p>' . nl2br($event['description']) . '</p>';
          } ?>
          <a href="<?php echo $event['url']; ?>" class="btn btn-outline-success">
            See details
          </a>
        <?php else : ?>
          <h6 class="small text-muted mb-0">
            <?php echo $date_string; ?> -
            <a class="text-success" href="<?php echo $event['url']; ?>">
              See details
            </a>
          </h6>
        <?php endif; ?>
      </div>
    </div>


  <?php
  endforeach;
}

//
// SINGLE EVENT
//
if (isset($_GET['event']) && substr($_GET['event'], 0, 7) == 'events/') {

  // Parse the markdown before header.php, so that we can override subtitle etc
  $markdown_fn = $md_base . $_GET['event'] . '.md';
  require_once('../includes/parse_md.php');

  //
  // Add event meta to the subtitle
  //

  $output = parse_md($markdown_fn);

  // Event type badge
  if (isset($meta['type'])) {
    $colour_class = $event_type_classes[strtolower($meta['type'])];
    $icon_class = $event_type_icons[strtolower($meta['type'])];
    $subtitle = '<span class="badge badge-' . $colour_class . ' mr-3"><i class="' . $icon_class . ' mr-1"></i>' . ucfirst($meta['type']) . '</span> ' . $subtitle;
  }

  $event = sanitise_date_meta($output["meta"]);

  if ($event) {
    $header_html = '<div class="row" style="margin-bottom:-1rem;"><div class="col-md-6">';
    $header_html .= '<dl>';
    // Start time
    if ($event['start_time']) {
      $header_html .= '<dt>Event starts:</dt><dd data-timestamp="' . $event['start_ts'] . '">' . date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) . '</dd>';
    } else {
      $header_html .= '<dt>Event starts:</dt><dd data-timestamp="' . $event['start_ts'] . '">' . date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) . '</dd>';
    }
    // End time
    if ($event['end_ts'] > $event['start_ts'] && $event['end_time']) {
      $header_html .= '<dt>Event ends:</dt><dd data-timestamp="' . $event['end_ts'] . '">' . date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']) . '</dd>';
    } else if ($event['end_ts'] > $event['start_ts']) {
      $header_html .= '<dt>Event ends:</dt><dd data-timestamp="' . $event['end_ts'] . '">' . date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']) . '</dd>';
    }
    $header_html .= '</dl>';
    $header_html .= '</div><div class="col-md-6">';
    // Location
    if (
      array_key_exists('location_name', $event) ||
      array_key_exists('location_url', $event) ||
      array_key_exists('address', $event) ||
      array_key_exists('location_latlng', $event)
    ) {
      if (isset($event['location_name'])) {
        $header_html .=  '<dt class="col-sm-3">Location:</dt><dd class="col-sm-9">';
        if (isset($event['location_url'])) {
          $header_html .=  '<a class="text-white underline" href="' . $event['location_url'] . '">' . $event['location_name'] . '</a>' . '<br>';
        } else {
          $header_html .=  $event['location_name'] . '<br>';
        }
      } else if (isset($event['location_url'])) {
        $header_html .=  '<dt>Web address:</dt><dd>';
        if (is_array($event['location_url'])) {
          foreach ($event['location_url'] as $url) {
            $header_html .=  '<a class="text-white underline" href="' . $url . '">' . $url . '</a>' . '<br>';
          }
        } else {
          $header_html .=  '<a class="text-white underline" href="' . $event['location_url'] . '">' . $event['location_url'] . '</a>' . '<br>';
        }
      }
      if (isset($event['address'])) {
        $header_html .=  $event['address'] . '<br>';
      }
      if (isset($event['location_latlng'])) {
        $header_html .=  '<a class="mt-2 btn btn-sm btn-outline-light" href="https://www.google.com/maps/search/?api=1&query=' . implode(',', $event['location_latlng']) . '" target="_blank">See map</a>';
      }
      $header_html .= '</dd>';
    }
    $header_html .= '</div></div>';
  }

  $md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/' . $_GET['event'] . '.md';

  // header.php runs parse_md() again to produce main page content
  $import_moment = true;
  $no_print_content = true;
  $mainpage_container = false;
  // Add in a YouTube embed if we have one
  if (array_key_exists('youtube_embed', $event)) {
    if (!is_array($event['youtube_embed'])) $event['youtube_embed'] = [$event['youtube_embed']];
    $youtube_embed = true;
    foreach ($event['youtube_embed'] as $embed) {
      $video_id = get_youtube_id($embed);
      if ($video_id) {
        echo'<script>var video_id="' . $video_id . '"</script>';
      }
    }
  }

  include('../includes/header.php');

  

  $toc = generate_toc($content);

  # only add ToC if there are more than two items in it
  if(substr_count($toc, "list-group-item ")>2){
    # Make a row with a column for content
    echo '<div class="row flex-wrap-reverse flex-lg-wrap"><div class="col-12 col-lg-9">';

    # Print content

    echo '<div class="rendered-markdown container container-xl main-content ml-5 pr-5">' . $content . '</div>';
    echo '</div>'; # close column div
    echo '<div class="col-12 col-lg-3 pl-2"><div class="side-sub-subnav sticky-top">';

    #add  ToC
    $toc = '<nav class="toc">' . $toc;

    # Add on the action buttons for the parameters docs

    # Back to top link
    $toc .= '<p class="small text-right"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
    $toc .= '</nav>';
    echo $toc;

    echo '</div></div>'; # end of the sidebar col
    echo '</div>'; # end of the row
  } else{
    echo '<div class="container main-content">';

    # Print content

    echo '<div class="rendered-markdown ">' . $content . '</div></div>';
  }

  // Javascript for moment time zone support
  if ($event['start_time']) {
    echo '
    <script type="text/javascript">
    $("[data-timestamp]").each(function(){
      var timestamp = $(this).data("timestamp");
      var timeformat = $(this).data("timeformat") ? $(this).data("timeformat") : "HH:mm z, LL";
      var local_time = moment.tz(timestamp, "X", moment.tz.guess());
      $(this).text(local_time.format(timeformat));
      if(moment(timestamp,"X").diff(moment().format())<(30*60*1000)){
        $(this).parent("tr").addClass("table-success"); // highlight row in schedule if current time is less than 30 minutes after time in row
      }
    });
    </script>
    ';
  }

  include('../includes/footer.php');
  exit;
}




//
// EVENTS LISTING PAGE
//

$title = 'Events';
$subtitle = 'Details of past and future nf-core meetups.';
$md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/events';
$header_btn_url = 'https://nf-co.re/events/rss';
$header_btn_text = '<i class="fas fa-rss mr-1"></i> RSS';

// Load event front-matter
$events = [];
$year_dirs = glob($md_base . 'events/*', GLOB_ONLYDIR);
foreach ($year_dirs as $year) {
  $event_mds = glob($year . '/*.md');
  foreach ($event_mds as $event_md) {
    // Load the file
    $md_full = file_get_contents($event_md);
    if ($md_full !== false) {
      $fm = parse_md_front_matter($md_full);
      // Add the URL
      $fm['meta']['url'] = '/events/' . basename($year) . '/' . str_replace('.md', '', basename($event_md));
      // Add to the events array
      $events[] = $fm['meta'];
    }
  }
}

# Parse dates and sort events by date
$future_events = [];
$past_events = [];
$current_events = [];
foreach ($events as $idx => $event) {

  $event = sanitise_date_meta($event);
  if (!$event) {
    unset($events[$idx]);
    continue;
  }

  # Update arrays
  if ($event['start_ts'] > time()) {
    $future_events[$idx] = $event;
  } else  if ($event['start_ts'] < time() && $event['end_ts'] > time()) {

    $current_events[$idx] = $event;
  } else {
    $past_events[$idx] = $event;
  }
}

# Sort future events so that the oldest is at the top
usort($future_events, function ($a, $b) {
  return $a['start_ts'] - $b['start_ts'];
});
# Sort past events so that the newest is at the top
usort($past_events, function ($a, $b) {
  return $b['start_ts'] - $a['start_ts'];
});

//
// RSS FEED
//
if (isset($_GET['rss'])) {
  header('Content-type: application/rss+xml');
  echo '<?xml version="1.0" encoding="UTF-8" ?>
    <rss version="2.0" xmlns:atom="http://www.w3.org/2005/Atom">
    <channel>
      <title>nf-core: ' . $title . '</title>
      <link>https://www.nf-co.re/events</link>
      <atom:link href="https://www.nf-co.re/events/rss" rel="self" type="application/rss+xml" />
      <description>' . $subtitle . '</description>
      ';
  if (count($future_events) > 0) {
    foreach ($future_events as $event) {
      echo '
      <item>
        <title>' . htmlspecialchars(utf8_encode($event['title'])) . '</title>
        <link>https://nf-co.re' . $event['url'] . '</link>
        <guid>https://nf-co.re' . $event['url'] . '</guid>
        <pubDate>' . date('r', $event['start_ts']) . '</pubDate>
        <description>' . htmlspecialchars(utf8_encode($event['subtitle'])) . '</description>
      </item>
      ';
    }
  }
  echo "\n    </channel>\n</rss>";
  exit;
}

//
// Web listing page
//

include('../includes/header.php');
echo '<div class="event-list">';
if (count($current_events) > 0) {
  echo '<h2 id="current_events"><a href="#current_events" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a><i class="fad fa-calendar mr-2"></i> Ongoing Events</h2>';
  print_current_events($current_events, true);
  echo '<hr>';
}

echo '<h2 id="future_events"><a href="#future_events" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a><i class="fad fa-calendar-day mr-2"></i> Upcoming Events</h2>';
if (count($future_events) > 0) {
  print_events($future_events, false);
} else {
  print '<p class="text-muted">No events found</p>';
}

echo '<hr>';
echo '<h2 id="past_events"><a href="#past_events" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a><i class="fad fa-calendar-check mr-2"></i> Past Events</h2>';
if (count($past_events) > 0) {
  print_events($past_events, true);
} else {
  print '<p class="text-muted">No events found</p>';
}
echo '</div>';
include('../includes/footer.php');
