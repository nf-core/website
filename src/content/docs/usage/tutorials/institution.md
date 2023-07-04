---
title: How to add yourself to the community page
subtitle: Adding your institution to the contributor's list
---

In case you couldn't find your institution / group of contributors on the [community pages](https://nf-co.re/contributors), please add yourself. All you need to add is a few lines of YAML and your institution's logo.

1. Fork the repository [nf-core/nf-core](https://github.com/nf-core/nf-co.re) to your own GitHub account. Look [here](https://guides.github.com/activities/forking/) for advice on how to fork a GitHub project.

2. In your own fork, modify the YAML file `nf-core-contributors.yaml` in the text editor of your choice.

3. Fill in the details about you and your institution directly into the YAML file:

   - `full_name:` Full Name of the Institution
   - `short_name:` Short Name of the Institution (e.g. an acronym)
   - `description:` Short description of the main activities of the Institution (see [`nf-core-contributors.yaml`](https://github.com/nf-core/nf-co.re/blob/master/nf-core-contributors.yaml) for examples)
   - `address:` Postal Address
   - `url:` Institution URL
   - `affiliation:` Affiliation - could be the subunit or department of your institution where you work (if _Fancy University_ is the Institution, then Affiliation would be the _McFancyDepartment_)
   - `affiliation_url:` Affiliation URL
   - `image_fn:` Logo filename (see below)
   - `contact`: Your Name
   - `contact_email:` Your E-Mail
   - `contact_github:` Your GitHub username
   - `location:` Institute location in geocoordinates format (`[<longitude>, <latitude>]`)
   - `twitter:` URL of the institutional twitter account, if available

4. Next, add two versions of your institutional logo in SVG format with the name that you have described in the YAML file under `image_fn:`

   - Upload a _white_ version to the folder `public_html/assets/img/contributors-white/`. It must be a single monochrome shape with no background colour.
   - Upload a _colour_ version to the folder `public_html/assets/img/contributors-colour/`.
   - See the website repository for examples in [white](https://github.com/nf-core/nf-co.re/tree/master/public_html/assets/img/contributors-white) and [colour](https://github.com/nf-core/nf-co.re/tree/master/public_html/assets/img/contributors-colour).
   - Both images should have the **same name** - if you define `image_fn: foobar.svg` in the YAML file, then both files should be named `foobar.svg` in the respective folders.
   - If you have only raster images available, please search for a SVG version. Ask for help in the nf-core Slack if in doubt (please do not add SVG files with embedded raster images - vector only).
     - _Tip_: Wikipedia often uses the SVG format of displayed logos available for the download.
     - If you were not successful at this point, skip it and let us know about it in the pull request.
     - _Note:_ Please make sure that you have the permission to use logos of your institution in an open source project. Organizations often don't have any concerns about it, and just want to be notified.

5. After you did the hardest part, please open a pull request from your fork. The modifications shall be compared to the master branch [nf-core/nf-core](https://github.com/nf-core/nf-co.re). We will provide you feedback if anything looks weird or hasn't been properly done. Feel free to let us know on Slack in the [`#request-review` channel](https://nfcore.slack.com/archives/CQY2U5QU9).
