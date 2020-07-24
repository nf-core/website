<?php
// Common PHP functions for the website

// Pull out YAML front-matter from a markdown file
require_once(dirname(__FILE__).'/libraries/Spyc.php');
function parse_md_front_matter($md_full){
    if(substr($md_full,0,3) == '---'){
      $md_parts = explode('---', $md_full, 3);
      if(count($md_parts) == 3){
        $meta = spyc_load($md_parts[1]);
        $md = $md_parts[2];
        return array(
            'meta' => $meta,
            'md' => $md
        );
      }
    }
    return array(
        'meta' => null,
        'md' => $md_full
    );
}

// From https://stackoverflow.com/a/18891474/713980
function time_ago($date, $ago=true) {
    $periods = array("second", "minute", "hour", "day", "week", "month", "year", "decade");
    $lengths = array("60", "60", "24", "7", "4.35", "12", "10");
    $now = time();
    if(is_numeric($date)) $unix_date = $date;
    else $unix_date = strtotime($date);
    // check validity of date
    if (empty($unix_date)) {
        return $date;
    }
    // is it future date or past date
    if ($now > $unix_date) {
        $difference = $now - $unix_date;
        $tense = "ago";
    } else {
        $difference = $unix_date - $now;
        $tense = "from now";
    }
    for ($j = 0; $difference >= $lengths[$j] && $j < count($lengths) - 1; $j++) {
        $difference /= $lengths[$j];
    }
    $difference = round($difference);
    if ($difference != 1) {
        $periods[$j].= "s";
    }
    $returnstring = "$difference $periods[$j]";
    if($ago || (!$ago && $tense != 'ago')){
        $returnstring .= " {$tense}";
    }
    return $returnstring;
}


function rsort_releases($a, $b){
    $t1 = strtotime($a->published_at);
    $t2 = strtotime($b->published_at);
    return $t2 - $t1;
}
function rsort_pipelines($a, $b){
    $t1 = strtotime($a->last_release);
    $t2 = strtotime($b->last_release);
    return $t2 - $t1;
}

function round_nicely($num){
  if($num > 1000000){
    $num /= 1000000;
    $num = round($num, 2).'M';
  } else if($num > 1000){
    $num /= 1000;
    $num = round($num, 2).'K';
  }
  return $num;
}

function endswith($haystack, $needle){
  $length = strlen( $needle );
  if( !$length ) { return true; }
  return substr( $haystack, -$length ) === $needle;
}

function return_json($response){
    // Spit out a JSON response with correct headers and exit
    header('Content-type: application/json');
    echo json_encode($response, JSON_PRETTY_PRINT);
    exit;
}

function get_self_url($strip_query=true){
    // Build URL for this page
    if(isset($_SERVER['HTTPS']) && $_SERVER['HTTPS'] === 'on') $self_url = "https://";
    else $self_url = "http://";
    if($strip_query){
      $url = strtok($_SERVER["REQUEST_URI"], '?');
    } else {
      $url = $_SERVER["REQUEST_URI"];
    }
    return $self_url.$_SERVER['HTTP_HOST'].$url;
}

function generate_toc($html_string){
  $toc = '';
  $curr_level = 0;
  $id_regex = "~<h([1-3])([^>]*)id\s*=\s*['\"]([^'\"]*)['\"]([^>]*)>(.*)</h[1-3]>~Uis";
  preg_match_all($id_regex, $html_string, $matches, PREG_SET_ORDER);
  if($matches){
    foreach($matches as $match){
      $whole_str = $match[0];
      $level = $match[1];
      $before_attrs = trim($match[2]);
      $id = trim($match[3]);
      $after_attrs = trim($match[4]);
      $h_content = $match[5];
      $name = trim(str_replace('&nbsp;','', htmlentities(strip_tags($h_content) )));
      if($level > $curr_level){
        $toc .= "\n".'<div class="list-group">'."\n";
      } else if($level == $curr_level) {
        $toc .= "\n";
      } else {
        while($level < $curr_level){
          $toc .= "\n</div>\n\n";
          $curr_level -= 1;
        }
      }
      $curr_level = $level;
      if(preg_match('/<code>.*?<\/code>/',$whole_str)){
        $name = '<code>'.$name.'</code>';
      }
      if(preg_match('/<i.*?<\/i>/',$whole_str,$icon_match)){
        $name = $icon_match[0].$name;
      }
      $name = str_replace('hidden','',$name); // remove artifact from "hidden" badge
      $is_hidden = strpos($before_attrs, 'toc-hidden') !== false || strpos($after_attrs, 'toc-hidden') !== false;
      $toc_hidden = $is_hidden ? 'collapse' : '';
      $toc .= '<a class="list-group-item list-group-item-action '.$toc_hidden.'" href="#'.$id.'">'.$name.'</a>';
    }
  }
  while($curr_level > 0){
    $toc .= '</div>';
    $curr_level -= 1;
  }
  return $toc;
}

$heading_ids = [];
function _h($level, $html, $toc_hidden=false){
  ////////////////
  // Build a heading tag with ID and anchor link
  ////////////////
  global $heading_ids;
  # Clean up the ID
  $hid = trim(strip_tags($html));
  $hid = strtolower( preg_replace('/[^\w\-\.]/', '', str_replace(' ', '-', $hid)));
  # Avoid duplicate IDs
  $i = 1; $base_hid = $hid;
  while(in_array($hid, $heading_ids)){
    $hid = $base_hid.'-'.$i;
    $i += 1;
  }
  # Class for hiding in ToC
  $toc_hidden_class = $toc_hidden ? 'toc-hidden' : '';
  return '
    <h'.$level.' id="'.$hid.'" class="'.$toc_hidden_class.'">
      <a href="#'.$hid.'" class="header-link"><span class="fas fa-link"></span></a>
      '.$html.'
    </h'.$level.'>';
};
function _h1($html){ return _h(1, $html); }
function _h2($html){ return _h(2, $html); }
function _h3($html){ return _h(3, $html); }
function _h4($html){ return _h(4, $html); }
function _h5($html){ return _h(5, $html); }


function add_ids_to_headers($content_input){
  //////////////////
  // Add IDs and anchor links to all headings in a block of HTML
  //////////////////
  global $heading_ids;
  $content_output = preg_replace_callback(
    '~<h([1234])>(.*?)</h([1234])>~Ui', // Ungreedy by default, case insensitive
    function ($matches) use($heading_ids) {
      $id_match = strip_tags($matches[2]);
      $id_match = strtolower( preg_replace('/[^\w\-\.]/', '', str_replace(' ', '-', $id_match)));
      $id_match = str_replace('---', '-', $id_match);
      $hid = $id_match;
      $i = 1;
      while(in_array($hid, $heading_ids)){
        $hid = $id_match.'-'.$i;
        $i += 1;
      }
      $heading_ids[] = $hid;
      return '<h'.$matches[1].' id="'.$hid.'"><a href="#'.$hid.'" class="header-link"><span class="fas fa-link"></span></a>'.$matches[2].'</h'.$matches[3].'>';
    },
    $content_input
  );
  return $content_output;
}
