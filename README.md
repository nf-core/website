
<img src="public_html/assets/img/logo/nf-core-logo.png" width="400">

# [nf-co.re](https://github.com/nf-core/nf-co.re)

This repository contains code for the nf-core website: **http://nf-co.re/**

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
To make edits to the website, fork the repository to your own user on GitHub and then clone to your local system.

**IMPORTANT:** The repo has git submodules, so remember to use the `--recursive` flag:

```
git clone --recursive git@github.com:[USERNAME]/nf-co.re.git
```

If you forget the recursive flag (I always do), the markdown conversion won't work. You can pull the submodules when you realise this with the following command:

```
git submodule update --init --recursive
```

To run the website locally, you need a standard AMP stack: Apache, MySQL and PHP (MySQL not needed at time of writing). For this, I recommend using the free version of [MAMP](https://www.mamp.info/en/).

Set the base directory to `/path/to/nf-co.re/public_html` in _Preferences > Web-Server > Document Root_ and then hit _Start Servers_.

I've built the site so that most of the hand-written text is in `/markdown`, to make it easier to write. The PHP files in `/public_html` then parse this into HTML dynamically, if supplied with a filename.

Note that the `.htaccess` file is set up to remove the `.php` file extensions in URLs.

## Credits
Phil ([@ewels](http://github.com/ewels/)) built this site, mostly over the course of one caffeine-fuelled evening.

## Help
If you have any questions or issues please send us a message on [Slack](https://nf-core-invite.herokuapp.com/).
