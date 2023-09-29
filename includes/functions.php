<?php
// Common PHP functions for the website

// Include the PHP Composer autoloader for required libraries
require dirname(__FILE__) . '/../vendor/autoload.php';

// Pull out YAML front-matter from a markdown file
function parse_md_front_matter($md_full) {
    if (substr($md_full, 0, 3) == '---') {
        $md_parts = explode('---', $md_full, 3);
        if (count($md_parts) == 3) {
            $meta = spyc_load($md_parts[1]);
            $md = $md_parts[2];
            return [
                'meta' => $meta,
                'md' => $md,
            ];
        }
    }
    return [
        'meta' => null,
        'md' => $md_full,
    ];
}
function get_correct_text_color($bg_color) {
    $text_color = in_array($bg_color, ['warning', 'light']) ? 'text-dark' : '';
    return $text_color;
}

// Helper function for event page
function sanitise_date_meta($event) {
    # Check that start date is set, delete if not
    if (!isset($event['start_date'])) {
        return false;
    }
    # Check end date is set
    if (!isset($event['end_date'])) {
        $event['end_date'] = $event['start_date'];
    }
    # Parse dates
    if (!isset($event['start_time'])) {
        $event['start_time'] = '';
    }
    if (!isset($event['end_time'])) {
        $event['end_time'] = '';
    }
    $event['start_ts'] = strtotime($event['start_date'] . ' ' . $event['start_time']);
    $event['end_ts'] = strtotime($event['end_date'] . ' ' . $event['end_time']);
    # Check end is after start
    if ($event['end_ts'] < $event['start_ts']) {
        $event['end_date'] = $event['start_date'];
        $event['end_ts'] = strtotime($event['end_date'] . ' ' . $event['end_time']);
    }
    return $event;
}

function prep_current_event($event) {
    $d = [];
    $d['event_type_classes'] = [
        'hackathon' => 'primary',
        'talk' => 'success',
        'poster' => 'secondary',
        'tutorial' => 'info',
        'training' => 'warning',
    ];
    $d['event_type_icons'] = [
        'hackathon' => 'fad fa-laptop-code',
        'talk' => 'fad fa-presentation',
        'poster' => 'fad fa-image',
        'tutorial' => 'fad fa-graduation-cap',
        'training' => 'fad fa-chalkboard-teacher',
    ];
    # Nice date strings
    $d['date_string'] =
        date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']) . ' - ' . date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    if (date('mY', $event['start_ts']) == date('mY', $event['end_ts'])) {
        $d['date_string'] =
            date('H:i', $event['start_ts']) .
            '-' .
            date('H:i e', $event['end_ts']) .
            ', ' .
            date('j<\s\u\p>S</\s\u\p> ', $event['start_ts']) .
            ' - ' .
            date('j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    }
    if (date('dmY', $event['start_ts']) == date('dmY', $event['end_ts'])) {
        $d['date_string'] =
            date('H:i', $event['start_ts']) . '-' . date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['end_ts']);
    }
    # Nice date strings
    if ($event['start_time']) {
        $d['nice_date_string'] = [
            'data-timestamp="' . $event['start_ts'] . '"',
            date('H:i e, j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']),
        ];
    } else {
        $d['nice_date_string'] = [
            'data-timestamp="' . $event['start_ts'] . '"',
            date('j<\s\u\p>S</\s\u\p> M Y', $event['start_ts']),
        ];
    }
    $d['colour_class'] = $d['event_type_classes'][strtolower($event['type'])];
    $d['text_colour_class'] = get_correct_text_color($d['colour_class']);
    $d['icon_class'] = $d['event_type_icons'][strtolower($event['type'])];
    $d['event_type_badge'] =
        '<span class="badge bg-' .
        $d['colour_class'] .
        ' ' .
        $d['text_colour_class'] .
        ' small"><i class="' .
        $d['icon_class'] .
        ' me-1"></i>' .
        ucfirst($event['type']) .
        '</span>';
    $d['location_url_meta'] = [];
    if (array_key_exists('location_url', $event)) {
        if (!is_array($event['location_url'])) {
            $event['location_url'] = [$event['location_url']];
        }
        foreach ($event['location_url'] as $idx => $url) {
            $d['location_url_meta'][$idx]['base_url'] = substr($url, 8, 7);
            switch ($d['location_url_meta'][$idx]['base_url']) {
                case 'zoom.us':
                    $d['location_url_meta'][$idx]['icon'] = '<i class="fas fa-video me-1"></i>';
                    $d['location_url_meta'][$idx]['print_url'] = count($event['location_url']) > 3 ? '' : $url;
                    break;
                case 'youtu.b':
                case 'youtube':
                    $d['location_url_meta'][$idx]['icon'] = '<i class="fab fa-youtube me-1"></i>';
                    $d['location_url_meta'][$idx]['print_url'] = count($event['location_url']) > 3 ? '' : $url;
                    break;
                case 'www.bil':
                    $d['location_url_meta'][$idx]['icon'] = '<i class="fad fa-film me-1"></i>';
                    $d['location_url_meta'][$idx]['print_url'] = count($event['location_url']) > 3 ? '' : $url;
                    break;
                case 'doi.org':
                    $d['location_url_meta'][$idx]['icon'] = '<i class="fas fa-barcode me-1"></i>';
                    $d['location_url_meta'][$idx]['print_url'] = count($event['location_url']) > 3 ? '' : $url;
                    break;
                default:
                    $d['location_url_meta'][$idx]['icon'] = '<i class="fas fa-external-link me-1"></i>';
                    $d['location_url_meta'][$idx]['print_url'] = $url;
            }
        }
    }
    return $d;
}

function print_current_events($events, $border) {
    foreach ($events as $idx => $event):

        $d = prep_current_event($event);
        $header_html .=
            '<dt>Event starts:</dt><dd ' . $d['nice_date_string'][0] . '>' . $d['nice_date_string'][1] . '</dd>';
        ?>

    <!-- Event Card -->
    <div class="card mb-3 <?php echo $border
        ? 'border-top-0 border-end-0 border-bottom-0 rounded-0 border-' . $d['colour_class']
        : 'border-0'; ?> ">
      <div class="card-body py-3 d-flex">
        <div class="pt-2"><i class="<?php echo $d['icon_class']; ?> fa-5x text-<?php echo $d[
     'colour_class'
 ]; ?> me-2"></i></div>
        <div class="px-2 flex-grow-1 d-flex flex-column justify-content-between">
          <div>
            <h5 class=" my-0 py-0 d-flex">
              <a class="text-success flex-grow-1" href="<?php echo $event['url']; ?>"><?php echo $event['title']; ?></a>
              <small><?php echo $d['event_type_badge']; ?></small>
            </h5>
            <?php
            if (array_key_exists('subtitle', $event)) {
                echo '<span class="mb-0 text-mute"> ' . $event['subtitle'] . '</span>';
            }
            if (array_key_exists('description', $event)) {
                echo '<p>' . nl2br($event['description']) . '</p>';
            }
            ?>
          </div>
          <div class="d-md-flex justify-content-between align-items-end">
            <h6 class=""><?php echo $d['date_string']; ?></h6>

            <a href="<?php echo $event['url']; ?>" class="btn btn-outline-success">
              See details
            </a>
          </div>
        </div>
      </div>
      <?php if (array_key_exists('location_url', $event) && $event['location_url'][0] != '#') {
          echo '<div class="btn-group mt-1" role="group" aria-label="External links">';
          foreach ($event['location_url'] as $idx => $url) {
              $m = $d['location_url_meta'][$idx];
              echo '<a href="' .
                  $url .
                  '" class="btn btn-' .
                  $d['colour_class'] .
                  ' rounded-0" data-bs-toggle="tooltip" title="' .
                  $url .
                  '">' .
                  $m['icon'] .
                  $m['print_url'] .
                  '</a>';
          }
          echo '</div>';
      } ?>
    </div>

<?php
    endforeach;
}

// From https://stackoverflow.com/a/18891474/713980
function time_ago($date, $ago = true) {
    $periods = ['second', 'minute', 'hour', 'day', 'week', 'month', 'year', 'decade'];
    $lengths = ['60', '60', '24', '7', '4.35', '12', '10'];
    $now = time();
    if (is_numeric($date)) {
        $unix_date = $date;
    } else {
        $unix_date = strtotime($date);
    }
    // check validity of date
    if (empty($unix_date)) {
        return $date;
    }
    // is it future date or past date
    if ($now > $unix_date) {
        $difference = $now - $unix_date;
        $tense = 'ago';
    } else {
        $difference = $unix_date - $now;
        $tense = 'from now';
    }
    for ($j = 0; $difference >= $lengths[$j] && $j < count($lengths) - 1; $j++) {
        $difference /= $lengths[$j];
    }
    $difference = round($difference);
    if ($difference != 1) {
        $periods[$j] .= 's';
    }
    $returnstring = "$difference $periods[$j]";
    if ($ago || (!$ago && $tense != 'ago')) {
        $returnstring .= " {$tense}";
    }
    return $returnstring;
}

function rsort_releases($a, $b) {
    $t1 = strtotime($a->published_at);
    $t2 = strtotime($b->published_at);
    return $t2 - $t1;
}
function rsort_pipelines($a, $b) {
    $t1 = strtotime($a->last_release);
    $t2 = strtotime($b->last_release);
    return $t2 - $t1;
}

function round_nicely($num) {
    if ($num > 1000000) {
        $num /= 1000000;
        $num = round($num, 2) . 'M';
    } elseif ($num > 1000) {
        $num /= 1000;
        $num = round($num, 2) . 'K';
    }
    return $num;
}

function endswith($haystack, $needle) {
    $length = strlen($needle);
    if (!$length) {
        return true;
    }
    return substr($haystack, -$length) === $needle;
}

function return_json($response) {
    // Spit out a JSON response with correct headers and exit
    header('Content-type: application/json');
    echo json_encode($response, JSON_PRETTY_PRINT);
    exit();
}

function get_url_protocol() {
    if (((isset($_SERVER['HTTPS']) && $_SERVER['HTTPS'] === 'on')) 
     || endsWith($_SERVER['HTTP_HOST'], 'tol.sanger.ac.uk') ){
        $protocol = 'https://';
    }else{
        $protocol = 'http://';
    }
    return $protocol;
}

function get_self_url($strip_query = true) {
    // Build URL for this page
    $self_url = get_url_protocol();

    if ($strip_query) {
        $url = strtok($_SERVER['REQUEST_URI'], '?');
    } else {
        $url = $_SERVER['REQUEST_URI'];
    }
    return $self_url . $_SERVER['HTTP_HOST'] . $url;
}

function generate_toc($html_string) {
    $toc = '<div  class="d-none d-lg-block"><strong class="ms-3 d-inline-block w-100 text-secondary border-bottom">On this page</strong>
        <div  style="max-height: calc(100vh - 150px); overflow: auto;">';
    $toc_md = '<div class="dropdown d-block d-lg-none">
                <a href="#" class="btn btn-secondary-outline bg-body dropdown-toggle float-end border" data-bs-toggle="dropdown">On this page</a>
                <div class="dropdown-menu toc-md">';
    $is_active = true;
    $id_regex = "~<h([1-3])([^>]*)id\s*=\s*['\"]([^'\"]*)['\"]([^>]*)>(.*)</h[1-3]>~Uis";
    preg_match_all($id_regex, $html_string, $matches, PREG_SET_ORDER);
    if ($matches) {
        $counter = $curr_level = 0;
        $shift = min(array_column($matches, 1)) - 1; # get the highest heading level and shift levels to start from 1
        foreach ($matches as $match) {
            $whole_str = $match[0];
            $level = $match[1] - $shift;
            $before_attrs = trim($match[2]);
            $id = trim($match[3]);
            $after_attrs = trim($match[4]);
            $h_content = $match[5];
            $name = trim(
                str_replace(
                    ['&nbsp;', '&amp;'],
                    ['', '&'],
                    htmlentities(strip_tags($h_content, $allowed_tags = ['code'])),
                ),
            );
            if ($level > $curr_level) {
                $toc .= "\n" . '<nav class="nav flex-column flex-nowrap">' . "\n";
                $counter += 1;
            } elseif ($level == $curr_level) {
                $toc .= "\n";
            } else {
                while ($level < $counter) {
                    $toc .= "\n</nav>\n\n";
                    $counter -= 1;
                }
            }
            $curr_level = $level;
            if (preg_match('/<code>.*?<\/code>/', $whole_str, $code_match)) {
                $name = preg_replace('/--/', '&#8288;-&#8288;-&#8288;', $name);
                $name = html_entity_decode($name);
            }
            if (preg_match('/<i.*?<\/i>/', $whole_str, $icon_match)) {
                $name = $icon_match[0] . $name;
            }
            $is_hidden = strpos($before_attrs, 'toc-hidden') !== false || strpos($after_attrs, 'toc-hidden') !== false;
            $toc_hidden = $is_hidden ? ' collapse ' : '';
            $active = $is_active ? ' active ' : '';
            $is_active = false;
            if ($level <= 2) {
                $toc_md .=
                    '<a class="dropdown-item' . $active . $toc_hidden . '" href="#' . $id . '">' . $name . '</a>';
            }

            $toc .=
                '<a class="nav-link scroll_to_link py-1 ' .
                $toc_hidden .
                $active .
                '" href="#' .
                $id .
                '">' .
                $name .
                '</a>';
        }
    }
    while ($counter > 0) {
        $toc .= '</nav>';
        $counter -= 1;
    }
    $toc_md .= '<a><hr class="dropdown-divider"></a>';
    $toc_md .= '<!-- tock_md_button_placeholder -->';
    $toc_md .= '<a class="dropdown-item" href="#"><i class="fas fa-arrow-to-top"></i> Back to top</a>';
    $toc_md .= '</div></div>';
    $toc .= '</div></div>';
    $toc = $toc_md . $toc;
    return $toc;
}

$heading_ids = [];
function _h($level, $html, $toc_hidden = false) {
    ////////////////
    // Build a heading tag with ID and anchor link
    ////////////////
    global $heading_ids;
    # Clean up the ID
    $hid = trim(strip_tags($html));
    $hid = strtolower(preg_replace('/[^\w\-\.]/', '', str_replace(' ', '-', $hid)));
    # Avoid duplicate IDs
    $i = 1;
    $base_hid = $hid;
    while (in_array($hid, $heading_ids)) {
        $hid = $base_hid . '-' . $i;
        $i += 1;
    }
    # Class for hiding in ToC
    $toc_hidden_class = $toc_hidden ? 'toc-hidden' : '';
    return '
    <h' .
        $level .
        ' id="' .
        $hid .
        '" class="' .
        $toc_hidden_class .
        '">
      ' .
        $html .
        '
      <a href="#' .
        $hid .
        '" class="header-link"><span class="fas fa-link"></span></a>
    </h' .
        $level .
        '>';
}
function _h1($html) {
    return _h(1, $html);
}
function _h2($html) {
    return _h(2, $html);
}
function _h3($html) {
    return _h(3, $html);
}
function _h4($html) {
    return _h(4, $html);
}
function _h5($html) {
    return _h(5, $html);
}

function add_ids_to_headers($content_input, $is_hidden = false) {
    //////////////////
    // Add IDs and anchor links to all headings in a block of HTML
    //////////////////
    global $heading_ids;
    $content_output = preg_replace_callback(
        '~<h([1234])>(.*?)</h([1234])>~Ui', // Ungreedy by default, case insensitive
        function ($matches) use ($heading_ids, $is_hidden) {
            $id_match = trim(strip_tags($matches[2]));
            $id_match = strtolower(preg_replace('/[^\w\-\.]+/', '', str_replace(' ', '-', $id_match)));
            $id_match = str_replace('.', '-', $id_match); // periods break the js code, because they are not valid in selector ids
            $hid = $id_match;
            $i = 1;
            while (in_array($hid, $heading_ids)) {
                $hid = $id_match . '-' . $i;
                $i += 1;
            }
            $hid = preg_replace('/^[\s\-]+/', '', $hid); // remove dashes from start of string (e.g. for parameter)
            $heading_ids[] = $hid;
            $hidden_class = $is_hidden ? 'toc-hidden' : '';
            $to_return = '<h' . $matches[1] . ' id="' . $hid . '" class="' . $hidden_class . '">';
            if (strpos($matches[2], '<code>')) {
                $to_return .=
                    $to_return .
                    '<a href="#' .
                    $hid .
                    '" class="header-link parameter-link scroll_to_link"><span class="fas fa-link fa-xs me-2"></span></a>' .
                    $matches[2] .
                    '</h' .
                    $matches[3] .
                    '>';
            } else {
                $to_return .=
                    $to_return .
                    $matches[2] .
                    ' <a href="#' .
                    $hid .
                    '" class="header-link scroll_to_link"><span class="fas fa-link fa-xs ms-1"></span></a></h' .
                    $matches[3] .
                    '>';
            }
            return $to_return;
        },
        $content_input,
    );
    return $content_output;
}

function get_youtube_id($url) {
    // https://stackoverflow.com/questions/3392993/php-regex-to-get-youtube-video-id#comment11552053_6121972
    preg_match(
        "#(?<=v=)[a-zA-Z0-9-]+(?=&)|(?<=v\/)[^&\n]+(?=\?)|(?<=embed/)[^&\n]+|(?<=v=)[^&\n]+|(?<=youtu.be/)[^&\n]+#",
        $url,
        $matches,
    );
    if ($matches) {
        return $matches[0];
    }
    return false;
}

// Load event front-matter
$md_base = dirname(dirname(__FILE__)) . '/markdown/';
$events = [];
$year_dirs = glob($md_base . 'events/*', GLOB_ONLYDIR);
foreach ($year_dirs as $year) {
    // Markdown files
    $event_mds = glob($year . '/*.md');
    // Event subdirectories
    $event_dirs = glob($year . '/*', GLOB_ONLYDIR);
    foreach ($event_dirs as $event_dir) {
        if (is_file($event_dir . '/index.md')) {
            $event_mds[] = $event_dir . '/index.md';
        }
    }

    foreach ($event_mds as $event_md) {
        // Load the file
        $md_full = file_get_contents($event_md);
        if ($md_full !== false) {
            $fm = parse_md_front_matter($md_full);
            // Add the URL
            if (basename($event_md) == 'index.md') {
                $fm['meta']['url'] = '/events/' . basename($year) . '/' . basename(dirname($event_md));
            } else {
                $fm['meta']['url'] = '/events/' . basename($year) . '/' . str_replace('.md', '', basename($event_md));
            }
            // Add to the events array
            $events[] = $fm['meta'];
        }
    }
}

# Look to see if we have an upcoming / ongoing event to show and pick one
$curr_event = false;
$additional_ongoing = 0;
$additional_upcoming = 0;
foreach ($events as $idx => $event) {
    $event = sanitise_date_meta($event);
    if (!$event) {
        unset($events[$idx]);
        continue;
    }
    if (isset($event['start_announcement'])) {
        $time_window = $event['start_ts'] - strtotime($event['start_announcement']);
    } elseif ($event['end_ts'] - $event['start_ts'] > 3600 * 5) {
        $time_window = 86400 * 28; // show announcement 7 days ahead for full day events
    } else {
        $time_window = 86400;
    }
    if ($event['start_ts'] < time() + $time_window && $event['end_ts'] > time()) {
        $current_events[$idx] = $event;

        // Ongoing event
        if ($event['start_ts'] < time() && $event['end_ts'] > time()) {
            $event['ongoing'] = true;
            if (!$curr_event) {
                $curr_event = $event;
            }
            // If multiple events run now, take the one with most recent start time
            elseif ($event['start_ts'] < $curr_event['start_ts']) {
                $curr_event = $event;
            } else {
                $additional_ongoing++;
            }
        }
        // Upcoming event
        else {
            $event['ongoing'] = false;
            if (!$curr_event) {
                $curr_event = $event;
            }
            // If multiple events come up, take the one with earliest start time
            elseif ($event['start_ts'] < $curr_event['start_ts']) {
                $curr_event = $event;
            } else {
                $additional_upcoming++;
            }
        }
    }
}

// get all modules from database
function get_modules() {
    $config = parse_ini_file('../config.ini');
    $conn = mysqli_connect(
        $config['host'],
        $config['username'],
        $config['password'],
        $config['dbname'],
        $config['port'],
    );
    $sql = 'SELECT * FROM nfcore_modules ORDER BY LOWER(name)';
    $modules = [];
    if ($result = mysqli_query($conn, $sql)) {
        if (mysqli_num_rows($result) > 0) {
            while ($row = mysqli_fetch_array($result)) {
                $row['keywords'] = explode(';', $row['keywords']);
                $row['authors'] = explode(';', $row['authors']);
                $row['tools'] = json_decode($row['tools'], true);
                $row['input'] = json_decode($row['input'], true);
                $row['output'] = json_decode($row['output'], true);
                $modules[] = $row;
            }
            // Free result set
            mysqli_free_result($result);
        } else {
            echo 'Oops! Something went wrong. Please try again later.';
        }
    }
    return $modules;
}
