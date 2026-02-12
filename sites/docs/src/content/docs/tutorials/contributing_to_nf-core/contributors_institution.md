---
title: Add your institution to the contributor's list
subtitle: A guide for adding your institution to the contributor's list
shortTitle: Add your institution
---

nf-core is by design a collaborative effort and would not exist if its many dedicated contributors. nf-core contributors are welcome to add their institution to the [contributors page](https://nf-co.re/contributors) on the nf-core website.

To add your institution to the contributors page:

1. Fork the [nf-core/website](https://github.com/nf-core/website) GitHub repository.

   :::tip
   See [Contributing to a project](https://guides.github.com/activities/forking/) to learn how to contribute to a project through forking.
   :::

2. Add your details and your institutional details to [`contributors.yaml`](https://github.com/nf-core/website/blob/main/sites/main-site/src/config/contributors.yaml) on your fork:
   - `full_name:` Full name of your institution
   - `short_name:` Short name of your institution (e.g., an acronym)
   - `description:` Short description of the main activities of your institution (see [`contributors.yaml`](https://github.com/nf-core/website/blob/main/sites/main-site/src/config/contributors.yaml) for examples)
   - `address:` Postal address of your institution
   - `url:` Institution URL
   - `affiliation:` Your affiliation (i.e., the department of your institution where you work)
   - `affiliation_url:` Your affiliation URL
   - `image_fn:` Logo filename (see below)
   - `contact`: Your name
   - `contact_email:` Your e-mail
   - `contact_github:` Your GitHub username
   - `location:` Institute location in [geocoordinates format](https://support.google.com/maps/answer/18539?hl=en&co=GENIE.Platform%3DDesktop) (`[<longitude>, <latitude>]`)
   - `mastodon:` Mastodon handle (if available)
   - `bluesky:` Bluesky handle (if available)

3. Add a colour and a white SVG of your institutional logo to your fork:
   - Add a _white_ SVG with no background to the `/public/images/contributors/white/` folder.
   - Add a _colour_ SVG with no background to the `/public/images/contributors/colour/` folder.

     :::note
     Both SVG files must have the name that you defined in the `contributors.yaml` under `image_fn:` and must be vectors. Do not add a SVG file with an embedded raster image.
     :::

     :::note
     You are responsible for making sure that you have the permission to use your institutional logos.
     :::

     :::tip
     See the website repository for examples of [white](https://github.com/nf-core/website/tree/main/public/images/contributors/white) and [colour](https://github.com/nf-core/website/tree/main/public/images/contributors/colour) SVGs.
     :::

4. Open a pull request from your fork to the main branch of the [nf-core/website](https://github.com/nf-core/website).

5. Share your PR on the nf-core [`#request-review`](https://nfcore.slack.com/archives/CQY2U5QU9) Slack channel and wait for approval.
