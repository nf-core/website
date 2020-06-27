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
