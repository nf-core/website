<?php
// Let's make a logo!

$textstring = trim(str_replace('logo/', '', $_GET['t']));
if(!isset($_GET['f'])){
    $textstring = strtolower($textstring);
    $textstring = preg_replace("/[^a-z]/", '', $textstring);
}
$filename = 'nf-core-'.preg_replace("/[^a-z]/", '', $textstring).'.png';

if(strlen($textstring) == 0){
    header('HTTP/1.1 404 Not Found');
    include('404.php');
    die();
}

// Load the base image
$image = imagecreatefrompng("assets/img/logo/nf-core-repologo-base.png");

// Keep PNG transparency
imageAlphaBlending($image, true);
imageSaveAlpha($image, true);

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
$cropped_image = imagecrop($image, ['x' => 0, 'y' => 0, 'width' => max($text_width, $min_width), 'height' => 1000]);
imageAlphaBlending($cropped_image, true);
imageSaveAlpha($cropped_image, true);

// Make and destroy image
header("Content-type: image/png");
header('Content-Disposition: filename="nf-core-'.$filename.'.png"');
imagepng($cropped_image);
imagedestroy($cropped_image);
imagedestroy($image);
