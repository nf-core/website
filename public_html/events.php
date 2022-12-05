<?php

# To get parse_md_front_matter() and sanitise_date_meta() functions
require_once '../includes/functions.php';
use Spatie\CalendarLinks\Link;

$md_base = dirname(dirname(__FILE__)) . '/markdown/';
$event_type_classes = [
    'bytesize' => 'success',
    'hackathon' => 'primary',
    'poster' => 'danger',
    'talk' => 'success',
    'tutorial' => 'info',
    'training' => 'warning',
];
$event_type_icons = [
    'bytesize' => 'fas fa-apple-core',
    'hackathon' => 'fas fa-laptop-code',
    'poster' => 'far fa-image',
    'talk' => 'fas fa-presentation',
    'tutorial' => 'fas fa-graduation-cap',
    'training' => 'fas fa-chalkboard-teacher',
];

function create_event_download_button($event, $button_style) {
    $start = DateTime::createFromFormat('U', $event['start_ts']);
    $start->setTimezone(new DateTimeZone('Europe/Amsterdam'));
    $end = DateTime::createFromFormat('U', $event['end_ts']);
    $end->setTimezone(new DateTimeZone('Europe/Amsterdam'));
    $address = $event['address'] ? $event['address'] : '';
    $address = $event['location_url'] ? $event['location_url'] : $address; # prefer url over address
    $address = is_array($address) ? $address[0] : $address; # if multiple location urls are given, take the first one
    $link = Link::create($event['title'], $start, $end)
        ->description($event['subtitle'] ? $event['subtitle'] : '')
        ->address($address);

    $event_download_button =
        '<div class="dropup btn-group" role="group">
          <button type="button" class="btn ' .
        $button_style .
        ' dropdown-toggle" href="#" role="button" data-bs-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
            <i class="far fa-calendar-plus me-1"></i> Export event
          </button>
          <div class="dropdown-menu">
            <a class="dropdown-item" href="' .
        $link->ics() .
        '" target="_blank"> Download iCal Event</a>
            <a class="dropdown-item" href="' .
        $link->google() .
        '" target="_blank"> Add to Google Calendar</a>
            <a class="dropdown-item" href="' .
        $link->webOutlook() .
        '" target="_blank"> Add to Microsoft Outlook</a>
          </div>
        </div>';
    return $event_download_button;
}

function print_events($events, $is_past_event) {
    global $event_type_classes;
    global $event_type_icons;
    $current_year = date('Y');
    foreach ($events as $idx => $event):

        # Nice date strings
        $date_string =
            date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) .
            ' - ' .
            date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
        if (date('mY', $event['start_ts']) == date('mY', $event['end_ts'])) {
            $date_string =
                date('j<\s\u\p>S</\s\u\p> ', $event['start_ts']) .
                ' - ' .
                date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
        }
        if (date('dmY', $event['start_ts']) == date('dmY', $event['end_ts'])) {
            $date_string = date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
        }
        # if event title starts with bytesize change event type
        if (strpos($event['title'], 'Bytesize') === 0) {
            $event['type'] = 'bytesize';
        }

        if (($current_year != date('Y', $event['start_ts'])) & $is_past_event) {
            $current_year = date('Y', $event['start_ts']);
            echo _h3($current_year);
        }
        $colour_class = $event_type_classes[strtolower($event['type'])];
        $text_colour_class = get_correct_text_color($colour_class);
        $icon_class = $event_type_icons[strtolower($event['type'])];
        ?>

    <!-- Event Card -->
    <div class="card my-4 border-3 border-top-0 border-end-0 border-bottom-0 rounded-0 border-<?php echo $colour_class; ?> overflow-visible <?php echo $event[
     'type'
 ]; ?>">
      <div class="card-body <?php if ($is_past_event) {
          echo 'py-2';
      } ?>">
        <h5 class="my-0 py-0">
          <a class="text-success text-decoration-none" href="<?php echo $event['url']; ?>"><?php echo $event[
    'title'
]; ?></a>
          <small><span class="badge bg-<?php echo $colour_class .
              ' ' .
              $text_colour_class; ?> float-end small"><i class="<?php echo $icon_class; ?> me-1"></i><?php echo ucfirst(
     $event['type'],
 ); ?></span></small>
        </h5>
        <?php
        if (array_key_exists('subtitle', $event)) {
            $tm = $is_past_event ? 'text-muted' : '';
            echo '<p class="mb-0 ' . $tm . '">' . $event['subtitle'] . '</p>';
        }

        if (!$is_past_event): ?>
          <h6 class="small text-muted"><?php echo $date_string; ?></h6>

          <?php if (array_key_exists('description', $event)) {
              echo '<p>' . nl2br($event['description']) . '</p>';
          } ?>
          <div class="btn-group" role="group">
            <a href="<?php echo $event['url']; ?>" class="btn btn-outline-success">
              See details
            </a>
            <?php echo create_event_download_button($event, 'btn-outline-success'); ?>
          </div>
        <?php else: ?>
          <h6 class="small text-muted mb-0">
            <?php echo $date_string; ?> -
            <a class="text-success" href="<?php echo $event['url']; ?>">
              See details
            </a>
          </h6>
        <?php endif;
        ?>
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
    $markdown_fn = $md_base . $_GET['event'];
    if (is_file($markdown_fn . '.md')) {
        // Regular single-page event
        $markdown_fn .= '.md';
    } elseif (is_dir($markdown_fn) && file_exists($markdown_fn . '/index.md')) {
        // Nested event index page
        $href_url_prepend = basename($markdown_fn) . '/';
        $markdown_fn = $markdown_fn . '/index.md';
    }

    require_once '../includes/parse_md.php';

    $output = parse_md($markdown_fn);
    $event = sanitise_date_meta($output['meta']);
    if ($event) {
        $header_html = '<div class="row" style="margin-bottom:-1rem;"><div class="col-md-6">';
        $header_html .= '<dl>';
        // Start time
        if ($event['start_time']) {
            $header_html .=
                '<dt>Event starts:</dt><dd data-timestamp="' .
                $event['start_ts'] .
                '">' .
                date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) .
                '</dd>';
        } else {
            $header_html .=
                '<dt>Event starts:</dt><dd data-timestamp="' .
                $event['start_ts'] .
                '">' .
                date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) .
                '</dd>';
        }
        // End time
        if ($event['end_ts'] > $event['start_ts'] && $event['end_time']) {
            $header_html .=
                '<dt>Event ends:</dt><dd data-timestamp="' .
                $event['end_ts'] .
                '">' .
                date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']) .
                '</dd>';
        } elseif ($event['end_ts'] > $event['start_ts']) {
            $header_html .=
                '<dt>Event ends:</dt><dd data-timestamp="' .
                $event['end_ts'] .
                '">' .
                date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']) .
                '</dd>';
        }
        $header_html .= '</dl>';
        $header_html .= $event['end_ts'] > time() ? create_event_download_button($event, 'btn-outline-light') : '';
        $header_html .= '</div><div class="col-md-6">';
        // Location
        if (
            array_key_exists('location_name', $event) ||
            array_key_exists('location_url', $event) ||
            array_key_exists('address', $event) ||
            array_key_exists('location_latlng', $event)
        ) {
            if (isset($event['location_name'])) {
                $header_html .= '<dt class="col-sm-3">Location:</dt><dd class="col-sm-9">';
                if (isset($event['location_url'])) {
                    if (is_array($event['location_url'])) {
                        foreach ($event['location_url'] as $url) {
                            $location = count($event['location_url']) == 1 ? $event['location_name'] : $url;
                            $header_html .=
                                '<a class="text-white underline" href="' . $url . '">' . $location . '</a>' . '<br>';
                        }
                    } else {
                        $header_html .=
                            '<a class="text-white underline" href="' .
                            $event['location_url'] .
                            '">' .
                            $event['location_name'] .
                            '</a>' .
                            '<br>';
                    }
                } else {
                    $header_html .= $event['location_name'] . '<br>';
                }
            } elseif (isset($event['location_url'])) {
                $header_html .= '<dt>Web address:</dt><dd>';
                if (is_array($event['location_url'])) {
                    foreach ($event['location_url'] as $url) {
                        $header_html .= '<a class="text-white underline" href="' . $url . '">' . $url . '</a>' . '<br>';
                    }
                } else {
                    $header_html .=
                        '<a class="text-white underline" href="' .
                        $event['location_url'] .
                        '">' .
                        $event['location_url'] .
                        '</a>' .
                        '<br>';
                }
            }
            if (isset($event['address'])) {
                $header_html .= $event['address'] . '<br>';
            }
            if (isset($event['location_latlng'])) {
                $header_html .=
                    '<a class="mt-2 btn btn-sm btn-outline-light" href="https://www.google.com/maps/search/?api=1&query=' .
                    implode(',', $event['location_latlng']) .
                    '" target="_blank">See map</a>';
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
        if (!is_array($event['youtube_embed'])) {
            $event['youtube_embed'] = [$event['youtube_embed']];
        }
        $youtube_embed = true;
        foreach ($event['youtube_embed'] as $embed) {
            $video_id = get_youtube_id($embed);
            if ($video_id) {
                echo '<script>var video_id="' . $video_id . '"</script>';
            }
        }
    }

    // load the typeform javascript for embed forms, if we have one
    if (array_key_exists('import_typeform', $event)) {
        $import_typeform = true;
    }
    include '../includes/header.php';

    $toc = generate_toc($content);

    # only add ToC if there are more than two items in it
    if (substr_count($toc, 'list-group-item ') > 2) {
        # Make a row with a column for content
        echo '<div class="row flex-wrap-reverse flex-lg-wrap"><div class="col-12 col-lg-9">';

        # Print content

        echo '<div class="rendered-markdown container container-xl main-content ms-5 pe-5">' . $content . '</div>';
        echo '</div>'; # close column div
        echo '<div class="col-12 col-lg-3 ps-2 h-100 sticky-top"><div class="side-sub-subnav">';

        #add  ToC
        $toc = '<nav class="toc">' . $toc;

        # Add on the action buttons for the parameters docs

        # Back to top link
        $toc .=
            '<p class="small text-end mt-3 d-none d-lg-block"><a href="#" class="text-muted"><i class="fas fa-arrow-to-top"></i> Back to top</a></p>';
        $toc .= '</nav>';
        echo $toc;

        echo '</div></div>'; # end of the sidebar col
        echo '</div>'; # end of the row
    } else {
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
      var row_time = new Date(timestamp*1000);
      var now = new Date();
      if(row_time.getMonth() == now.getMonth() && row_time.getDate() <= now.getDate() && row_time.getHours()*60+row_time.getMinutes()<= now.getHours()*60+now.getMinutes() && now.getHours()*60+now.getMinutes() < row_time.getHours()*60+row_time.getMinutes()+30){
        $(this).parent("tr").addClass("table-success"); // highlight row in schedule if current time is less than 30 minutes after time in row
      }
    });
    </script>
    ';
    }

    include '../includes/footer.php';
    exit();
}

//
// EVENTS LISTING PAGE
//

$title = 'Events';
$subtitle = 'Details of past and future nf-core meetups.';
$md_github_url = 'https://github.com/nf-core/nf-co.re/tree/master/markdown/events';
$header_btn_url = 'https://nf-co.re/events/rss';
$header_btn_text = '<i class="fas fa-rss me-1"></i> RSS';

// Load event front-matter
$events = [];
$year_dirs = glob($md_base . 'events/*', GLOB_ONLYDIR);
foreach ($year_dirs as $year) {
    // Single page events
    $event_mds = glob($year . '/*');
    foreach ($event_mds as $fpath) {
        if (is_file($fpath) && substr($event_md, -3) == '.md') {
            $event_md = $fpath;
            $url = '/events/' . basename($year) . '/' . str_replace('.md', '', basename($event_md));
        } elseif (is_dir($fpath) && file_exists($fpath . '/index.md')) {
            $event_md = $fpath . '/index.md';
            $url = '/events/' . basename($year) . '/' . basename($fpath);
        }
        // Load the file
        $md_full = file_get_contents($event_md);
        if ($md_full !== false) {
            $fm = parse_md_front_matter($md_full);
            // Add the URL
            $fm['meta']['url'] = $url;
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
    } elseif ($event['start_ts'] < time() && $event['end_ts'] > time()) {
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
      <title>nf-core: ' .
        $title .
        '</title>
      <link>https://www.nf-co.re/events</link>
      <atom:link href="https://www.nf-co.re/events/rss" rel="self" type="application/rss+xml" />
      <description>' .
        $subtitle .
        '</description>
      ';
    if (count($future_events) > 0) {
        foreach ($future_events as $event) {
            echo '
      <item>
        <title>' .
                htmlspecialchars(utf8_encode($event['title'])) .
                '</title>
        <link>https://nf-co.re' .
                $event['url'] .
                '</link>
        <guid>https://nf-co.re' .
                $event['url'] .
                '</guid>
        <pubDate>' .
                date('r', $event['start_ts']) .
                '</pubDate>
        <description>' .
                htmlspecialchars(utf8_encode($event['subtitle'])) .
                '</description>
      </item>
      ';
        }
    }
    echo "\n    </channel>\n</rss>";
    exit();
}

//
// Web listing page
//

include '../includes/header.php';

// add a row of buttons
echo '<div class="btn-toolbar events-toolbar mb-4"><button type="button" class="btn txt-body">Filter:</button>';
echo '<div class="event-filters input-group input-group-sm me-2">
        <input type="search" class="form-control w-25" placeholder="Search events">
        </div>';
echo '<div class="btn-group events-filters w-50 align-items-center"  role="group" >';
echo '<input type="radio" class="btn-check" name="filter" id="filter-all" autocomplete="off" checked>
  <label class="btn btn-outline-secondary" for="filter-all">All</label>';

foreach ($event_type_classes as $class => $color) {
    echo '<input type="radio" class="btn-check " name="filter" data-bs-target=".' .
        $class .
        '" id="filter-' .
        $class .
        '" autocomplete="off">
  <label class="text-nowrap btn btn-outline-' .
        $color .
        '" for="filter-' .
        $class .
        '">' .
        '<i class=" ' .
        $event_type_icons[$class] .
        ' me-1"></i> ' .
        $class .
        '</label>';
}

echo '</div>';
echo '</div>';
echo '<div class="event-list">';
if (count($current_events) > 0) {
    echo '<div class="mb-5">';
    echo _h2('<i class="fad fa-calendar me-2"></i> Ongoing Events');
    print_current_events($current_events, true);
    echo '</div>';
    echo '<hr>';
}
echo '<div class="mb-5">';
echo _h2('<i class="fad fa-calendar-day me-2"></i> Upcoming Events');

if (count($future_events) > 0) {
    print_events($future_events, false);
} else {
    print '<p class="text-muted no-events">No events found</p>';
}
echo '</div>';
echo '<hr>';
echo '<div class="mb-5">';
echo _h2('<i class="fad fa-calendar-check me-2"></i> Past Events');
if (count($past_events) > 0) {
    print_events($past_events, true);
} else {
    print '<p class="text-muted no-events">No events found</p>';
}
echo '</div>';
echo '</div>';
echo '</div>';
include '../includes/footer.php';
