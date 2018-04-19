# nf-co.re

This repository contains code for the nf-core website: **http://nf-co.re**

<img src="public_html/assets/img/logo/nf-core-logo.png" width="500">

## Packages used
Here's how the website is built:

* Language: PHP
* HTML / CSS / JS framework: [Bootstrap v4](http://getbootstrap.com/)
* JavaScript libraries:
    * [jQuery](https://jquery.com/)
    * [Popper.js](https://popper.js.org/) _(used for bootstrap tooltips)_
    * [highlightjs](https://highlightjs.org/)
* PHP Markdown parsing: [Parsedown](https://github.com/erusev/parsedown/) and [Parsedown Extra](https://github.com/erusev/parsedown-extra/)
* SVG icons: http://www.flaticon.com, https://worldvectorlogo.com/

## Development
To run, the website needs a standard AMP stack: Apache, MySQL and PHP (MySQL not needed at time of writing).

For running locally, I recommend using [MAMP](https://www.mamp.info/en/) (OSX) / [WAMP](http://www.wampserver.com/en/) (Windows).

I've built the site so that most of the hand-written text is in `/markdown`, to make it easier to write. The PHP files in `/public_html` then parse this into HTML dynamically, if supplied with a filename.

Note that the `.htaccess` file is set up to remove the `.php` file extensions in URLs.

## Credits
Phil ([@ewels](http://github.com/ewels/)) built this site, mostly over the course of one caffeine-fuelled evening.
