<?php
// Let's make a logo!

$textstring = trim(str_replace('logo/', '', $_GET['t']));
if (!isset($_GET['f'])) {
    $textstring = strtolower($textstring);
    $textstring = preg_replace('/[^a-z]/', '', $textstring);
}

$new_width = false;
if (isset($_GET['w'])) {
    $new_width = $_GET['w'];
} elseif (isset($_GET['s'])) {
    $new_width = 400;
}

$theme = '';
if (isset($_GET['theme'])) {
    if ($_GET['theme'] == 'dark') {
        $theme = 'dark';
    } elseif ($_GET['theme'] == 'light') {
        $theme = 'light';
    }
}

$filename = 'nfcore-' . preg_replace('/[^a-z]/', '', $textstring) . '_logo_' . $theme . '.png';
$filename = str_replace('_.', '.', $filename);
$filename = basename($filename);

if (strlen($textstring) == 0) {
    header('HTTP/1.1 404 Not Found');
    include '404.php';
    die();
}

// Check if we have a cached version already
$cache_filename = $filename;
if ($new_width && is_numeric($new_width)) {
    $cache_filename =
        'nfcore-' . preg_replace('/[^a-z]/', '', $textstring) . '_logo_w' . $new_width . '_' . $theme . '.png';
    $cache_filename = str_replace('_.', '.', $cache_filename);
    $cache_filename = basename($cache_filename);
}
$logo_cache_fn = dirname(dirname(__FILE__)) . "/api_cache/logos/{$cache_filename}";
# Build directories if needed
if (!is_dir(dirname($logo_cache_fn))) {
    mkdir(dirname($logo_cache_fn), 0777, true);
}
// Return the cached version if it exists
if (file_exists($logo_cache_fn) && !isset($_GET['f'])) {
    header('Content-type: image/png');
    header('Content-Disposition: filename="' . $filename . '"');
    echo file_get_contents($logo_cache_fn);
    exit();
}

// Load the base image
if ($theme == 'dark') {
    $template_fn = 'assets/img/logo/nf-core-repo-logo-base-darkbg.png';
} else {
    $template_fn = 'assets/img/logo/nf-core-repo-logo-base-lightbg.png';
}
[$width, $height] = getimagesize($template_fn);
$image = imagecreatefrompng($template_fn);

// Create some colors
$black = imagecolorallocate($image, 0, 0, 0);
$color = $theme == 'dark' ? imagecolorallocate($image, 250, 250, 250) : imagecolorallocate($image, 5, 5, 5);
$font_size = 300;
$font_path = '../includes/Maven_Pro/MavenPro-Bold.ttf';

// Put text into image
imagettftext(
    $image, // image
    $font_size, // size
    0, // angle
    110, // x
    850, // y
    $color, // colour
    $font_path, // font
    $textstring, // text
);

// Crop off the excessive whitespace
$text_bbox = imagettfbbox($font_size, 0, $font_path, $textstring);
$text_width = abs($text_bbox[4] - $text_bbox[0]) + 250;
$min_width = 2300;
$width = max($text_width, $min_width);
$image = imagecrop($image, ['x' => 0, 'y' => 0, 'width' => $width, 'height' => $height]);

// If a width is given, scale the image
if (is_numeric($new_width)) {
    #$image = imagescale($image, 400, -1, IMG_NEAREST_NEIGHBOUR);
    // Get new dimensions
    $resize_factor = $new_width / $width;
    $new_height = $height * $resize_factor;

    // Create new image, with transparency
    $image_p = imagecreatetruecolor($new_width, $new_height);
    imagesetinterpolation($image_p, IMG_BICUBIC);
    imagealphablending($image_p, false);
    imagesavealpha($image_p, true);
    $transparent = imagecolorallocatealpha($image_p, 255, 255, 255, 127);

    // Resample
    imagecopyresampled($image_p, $image, 0, 0, 0, 0, $new_width, $new_height, $width, $height);

    // Overwrite full size image with resampled
    $image = $image_p;
}

// Keep PNG transparency
imageAlphaBlending($image, true);
imageSaveAlpha($image, true);

# Save image to cache
imagepng($image, $logo_cache_fn);
imagedestroy($image);

// Send the image to the browser
header('Content-type: image/png');
header('Content-Disposition: filename="' . $filename . '"');
echo file_get_contents($logo_cache_fn);

// Kill the cache if this was forced
if (isset($_GET['f'])) {
    unlink($logo_cache_fn);
}
