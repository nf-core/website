<img src="public_html/assets/img/logo/sanger-tol-logo.svg#gh-light-mode-only" width="400">
<img src="public_html/assets/img/logo/sanger-tol-logo-darkbg.svg#gh-dark-mode-only" width="400">

# [sanger-tol](https://pipelines.tol.sanger.ac.uk)

This repository contains code for the sanger-tol pipelines website: **<https://pipelines.tol.sanger.ac.uk>**. It has been forked and adapted from **<https://nf-co.re/>**.

## Packages used

Here's how the website is built:

- Language: PHP
- HTML / CSS / JS framework: [Bootstrap v5](http://getbootstrap.com/)
- JavaScript libraries:
  - [jQuery](https://jquery.com/)
  - [Popper.js](https://popper.js.org/) _(used for bootstrap tooltips)_
  - [highlightjs](https://highlightjs.org/) _(syntax highlighting)_
  - [Leaflet](https://leafletjs.com/) _(contributor map)_
  - [Moment.js](https://momentjs.com/) _(time and date parsing)_
  - [Chart.js](https://www.chartjs.org/) _(statistics plots)_
  - [hammer.js](https://hammerjs.github.io/) _(mobile touch interaction handling)_
  - [chartjs-plugin-zoom](https://github.com/chartjs/chartjs-plugin-zoom) _(Zoom and pan plugin for Chart.js)_
  - [Canvas2Svg.js](https://gliffy.github.io/canvas2svg/) _(SVG exports of Chart.JS plots)_
  - [FileSaver.js](https://github.com/eligrey/FileSaver.js/) _(Trigger browser downloads from in-page data, used to save plot SVGs to files)_
  - [jQuery table sorter](https://mottie.github.io/tablesorter/) _(sorting tables)_
- PHP Markdown parsing: [Parsedown](https://github.com/erusev/parsedown/) and [Parsedown Extra](https://github.com/erusev/parsedown-extra/)
- SVG icons: <http://www.flaticon.com>, <https://worldvectorlogo.com/>

## Development

### Getting the code

To make edits to the website, fork the repository to your own user on GitHub and then clone to your local system.

```bash
git clone git@github.com:[USERNAME]/pipelines-website.git
cd pipelines-website/
```

### Running a local server

Ok, you're ready! To run the website locally, just start the apache-php server with:

```bash
docker-compose up
```

> NB: If you are using a Mac with Apple silicon, you need to run:
>
> ```bash
> docker-compose -f m1-docker-compose.yml up
> ```

You should then be able to access the website in your browser at [http://localhost:8888/](http://localhost:8888/).

If you prefer, you can also use a tool such as [MAMP](https://www.mamp.info/) - if so,
set the base directory to `/path/to/pipelines-website/public_html` in _Preferences > Web-Server > Document Root_ and then hit _Start Servers_.

Most of the hand-written text is in `/markdown`, to make it easier to write. The PHP files in `/public_html` then parse this into HTML dynamically, if supplied with a filename.

Note that the `.htaccess` file is set up to remove the `.php` file extensions in URLs.

### First-run

Much of the site is powered by a `pipelines.json` file.
The webserver does this automatically when GitHub events trigger an update, but you'll need to run the script manually.

#### Access tokens

First you'll need a `config.ini` text file with values for `github_username` and `github_access_token`.
See [instructions on how to get a GitHub OAuth token](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line) (the token only needs the `public_repo` permission).
This file is ignored in `.gitignore` for security reasons.

For the MySQL database you should also add the following values:

```ini
host = 'db'
port = '3306';
dbname = 'nfcore';
username = 'nfcore_admin';
password = 'PEBBLY8exhibit_mead1cilium6despise'
```

#### Running PHP scripts

It's easiest to run these first manual update scripts on the command line. If you have PHP available
then you may be able to do this directly. Alternatively, if you are using Docker as above then you can
open a shell inside the running container. The container is typically named `nf-core-web` (you can check this
with the `docker ps` command), so you can open an interactive shell using the following command:

```bash
docker exec -it nf-core-web /bin/bash
cd /var/www/
```

#### Update scripts

The following command will create `public_html/pipelines.json`, which is used by the website.

```bash
php update_pipeline_details.php
```

To update the modules database (from within the docker container) run:

```bash
php update_module_details.php
```

Optionally, once you've done that, you can grab the pipeline traffic, issue statistics and font awesome icons:

```bash
php update_issue_stats.php
php update_stats.php
php update_fontawesome_icons.php
```

This creates `nfcore_stats.json`, `nfcore_issue_stats.json` and `public_html/assets/js/fa-icons.json`,
all also ignored in `.gitignore`.

## Production Server Setup

### Deployment

```bash
# login to the production VM

# clone the repo if not done yet
git clone https://github.com/sanger-tol/pipelines-website.git

# pull the latest from Git and the main branch being used.
cd pipelines-website
git status
git pull

# start the docker containers
docker-compose -f docker-compose-prod.yml up >> logs/docker-compose.log 2>&1 &
```

### Stats cronjob

The web server needs the following cronjobs running to scrape statistics and updates:

```cron
00 */2 * * * ( git -C /home/ubuntu/pipelines-website pull && docker run --rm -v /home/ubuntu/pipelines-website:/var/www/ -w /var/www composer:2.1 bash -c 'composer install' && docker run --rm --user 1000:1000 -v /home/ubuntu/pipelines-website:/var/www/ -w /var/www node:16-alpine3.12 sh -c 'npm install & npm run build-prod' ) >> /home/ubuntu/pipelines-website/logs/git-pull.log 2>&1
00 05 * * * docker exec nf-core-web /usr/local/bin/php /var/www/update_pipeline_details.php >> /home/ubuntu/pipelines-website/logs/update_pipelines.log 2>&1
15 05 * * * docker exec nf-core-web /usr/local/bin/php /var/www/update_module_details.php >> /home/ubuntu/pipelines-website/logs/update_modules.log 2>&1
30 05 * * * docker exec nf-core-web /usr/local/bin/php /var/www/update_stats.php >> /home/ubuntu/pipelines-website/logs/update_stats.log 2>&1
45 05 * * * docker exec nf-core-web /usr/local/bin/php /var/www/update_issue_stats.php >> /home/ubuntu/pipelines-website/logs/update_issue_stats.log 2>&1
```

Remember to replace the path with your actual deployment directory if any different.

The `update_issue_stats.php` script can use a lot of GitHub API calls, so should run at least one hour after the `update_stats.php` script last finished.
This is not because the script takes an hour to run, but because the GitHub API rate-limiting counts the number of calls within an hour.

## Contribution guidelines

If you are looking forward to contribute to the website or add your institution to the official list of contributors, please have a look at the [CONTRIBUTING.md](https://github.com/sanger-tol/pipelines-website/blob/main/.github/CONTRIBUTING.md).

## Community

If you have any questions or issues, please [let us know](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato%2Cmuffato&labels=connect&projects=&template=contact_us.yaml&title=%5BContact+Us%5D%3A+).

## Credits

Matthieu Muffato ([@muffato](http://github.com/muffato)) manages the content and Guoying Qi ([@gq1](https://github.com/gq1)) manages the website.
Many individuals, especially Priyanka Surana ([@priyanka-surana](http://github.com/priyanka-surana/)), have made various contributions.

Phil Ewels ([@ewels](http://github.com/ewels/)) built the original nf-core website.
More recently, [@mashehu](https://github.com/mashehu) has done a great deal of work with the code.
See the [repo contributors](https://github.com/nf-core/nf-co.re/graphs/contributors) for more.

Kudos to the excellent [npm website](https://www.npmjs.com), which provided inspiration for the design of the pipeline pages.
