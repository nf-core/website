<?php
// Let's make a logo!

$textstring = trim(str_replace('logo/', '', $_GET['t']));
if(!isset($_GET['f'])){
    $textstring = strtolower($textstring);
    $textstring = preg_replace("/[^a-z]/", '', $textstring);
}
$filename = 'nfcore-'.preg_replace("/[^a-z]/", '', $textstring).'_logo.png';

if(strlen($textstring) == 0){
    header('HTTP/1.1 404 Not Found');
    include('404.php');
    die();
}

// Load the base image
$image = imagecreatefrompng("assets/img/logo/nf-core-repologo-base.png");

// Create some colors
$black = imagecolorallocate($image, 0, 0, 0);
$font_size = 300;
$font_path = "../includes/Maven_Pro/MavenPro-Bold.ttf";

// Put text into image
imagettftext(
    $image,      // image
    $font_size,  // size
    0,           // angle
    110,         // x
    850,         // y
    $black,      // colour
    $font_path,  // font
    $textstring  // text
);

// Crop off the excessive whitespace
$text_bbox = imagettfbbox($font_size, 0, $font_path, $textstring);
$text_width = abs($text_bbox[4] - $text_bbox[0]) + 250;
$min_width = 2300;
$image = imagecrop($image, ['x' => 0, 'y' => 0, 'width' => max($text_width, $min_width), 'height' => 1000]);

// If a width is given, scale the image
if(isset($_GET['w']) && is_numeric($_GET['w'])){
    $image = imagescale($image, $_GET['w']);
} else if(isset($_GET['s'])){
    $image = imagescale($image, 400);
}

// Keep PNG transparency
imageAlphaBlending($image, true);
imageSaveAlpha($image, true);

// Make and destroy image
header("Content-type: image/png");
header('Content-Disposition: filename="'.$filename.'"');
imagepng($image);
imagedestroy($image);
imagedestroy($image);
