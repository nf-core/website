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
    * [highlightjs](https://highlightjs.org/) _(syntax highlighting)_
    * [Leaflet](https://leafletjs.com/) _(contributor map)_
    * [Moment.js](https://momentjs.com/) _(time and date parsing)_
    * [Chart.js](https://www.chartjs.org/) _(statistics plots)_
    * [hammer.js](https://hammerjs.github.io/) _(mobile touch interaction handling)_
    * [chartjs-plugin-zoom](https://github.com/chartjs/chartjs-plugin-zoom) _(Zoom and pan plugin for Chart.js)_
    * [Canvas2Svg.js](https://gliffy.github.io/canvas2svg/) _(SVG exports of Chart.JS plots)_
    * [FileSaver.js](https://github.com/eligrey/FileSaver.js/) _(Trigger browser downloads from in-page data, used to save plot SVGs to files)_
    * [jQuery table sorter](https://mottie.github.io/tablesorter/) _(sorting tables)_
* PHP Markdown parsing: [Parsedown](https://github.com/erusev/parsedown/) and [Parsedown Extra](https://github.com/erusev/parsedown-extra/)
* SVG icons: http://www.flaticon.com, https://worldvectorlogo.com/
## Development

To make edits to the website, fork the repository to your own user on GitHub and then clone to your local system.

**IMPORTANT:** The repo has git submodules, so remember to use the `--recursive` flag:

```bash
git clone --recursive git@github.com:[USERNAME]/nf-co.re.git
cd nf-co.re/
```

If you forget the recursive flag (I always do), the markdown conversion won't work. You can pull the submodules when you realise this with the following command:

```bash
git submodule update --init --recursive
```

Next, you'll need to build the `pipelines.json` file that powers much of the site. The webserver does this automatically when GitHub events trigger an update, but you'll need to run the script manually. Assuming you have PHP available on the command line, you can do this as follows.

First you'll need a `config.ini` text file with values for `github_username` and `github_access_token` set in it
([instructions on how to get a GitHub OAuth token](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line),
the token only needs the `public_repo` permission). This file is ignored in `.gitignore` for security reasons. Then you can run:

```bash
php update_pipeline_details.php
```

This will create `public_html/pipelines.json`, which is used by the website.
Note that this is ignored in the `.gitignore` file and will not be tracked in git history.

Optionally, once you've done that, you can grab the pipeline traffic and issue statistics:

```bash
php update_issue_stats.php
php update_stats.php
```

Note that your github account needs push rights for the nf-core permission for the `update_stats.php` to work.

This creates `nfcore_stats.json` and `nfcore_issue_stats.json`, also ignored in `.gitignore`.

Ok, you're ready! To run the website locally, you need a standard AMP stack: Apache, MySQL and PHP (MySQL not needed at time of writing). For this, I recommend using the free version of [MAMP](https://www.mamp.info/en/).

Set the base directory to `/path/to/nf-co.re/public_html` in _Preferences > Web-Server > Document Root_ and then hit _Start Servers_.

I've built the site so that most of the hand-written text is in `/markdown`, to make it easier to write. The PHP files in `/public_html` then parse this into HTML dynamically, if supplied with a filename.

Note that the `.htaccess` file is set up to remove the `.php` file extensions in URLs.

## Server Setup

### Stats cronjob

The web server needs the following cronjob running to scrape pipeline statistics once a week:

```
0	0	*	*	*	/usr/local/bin/php /home/nfcore/nf-co.re/update_stats.php >> /home/nfcore/update.log 2>&1
0	2	*	*	*	/usr/local/bin/php /home/nfcore/nf-co.re/update_issue_stats.php >> /home/nfcore/update.log 2>&1
```

The `update_issue_stats.php` script can use a lot of GitHub API calls, so should run at least one hour after the `update_stats.php` script last finished.

### Tools API docs

The repo has a softlink for `/tools-docs` which is intended for use on the server and corresponds to the path used in `public_html/deploy.php`. This script pulls the built API docs from the tools repo onto the server so that it can be served at that URL.

## Contribution guidelines

If you are looking forward to contribute to the website or add your institution to the official list of contributors, please have a look at the [CONTRIBUTING.md](./.github/CONTRIBUTING.md).

## Community

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## Credits

Phil Ewels ([@ewels](http://github.com/ewels/)) built the website, but there have been many contributors to the content and documentation. See the [repo contributors](https://github.com/nf-core/nf-co.re/graphs/contributors) for more.

Kudos to the excellent [npm website](https://www.npmjs.com), which provided inspiration for the design of the pipeline pages.
