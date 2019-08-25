<?php
// Common PHP functions for the website

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
