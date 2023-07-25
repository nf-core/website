---
title: Getting started
subtitle: How to run your first nf-core pipeline.
menu:
  main:
    weight: 10
---

# nf-core/website - Contributing Guidelines

Hi there! Many thanks for taking an interest in improving the nf-core website.

## Instructions to add your institution to the list of contributors

1. Fork this repo
2. Add your institute details and your name to the [contributors.yaml](../src/config/contributors.yaml)

Here's the desired format

```yaml
- full_name: Centre for Genomic Regulation
  short_name: CRG
  description: >
    The Centre for Genomic Regulation (CRG) is an international biomedical research institute
    of excellence, based in Barcelona, Spain, whose mission is to discover and advance knowledge
    for the benefit of society and public health.
  address: Carrer Dr. Aiguader, 88, 08003 Barcelona, Spain
  url: http://www.crg.eu/
  image_fn: CRG.svg
  contact: Paolo Di Tommaso
  contact_email: paolo.ditommaso@crg.eu
  contact_github: pditommaso
  location: [41.3853828, 2.191863]
  twitter: CRGenomica
```

3. For the purpose of displaying your institution logo, please also add the `SVG` of your institutions logo. These logos are displayed on the [main page](https://nf-co.re/) as
   well as in the [community page](https://nf-co.re/contributors#organisations).

Please note that`nf-co.re` website makes use of the `light/dark` mode so it's recommended that you add the logo for both modes. You could make use of `Adobe Illustrator/Inkscape` for converting the logo in two variants.

- The `colored SVG logo` needs to be added to the [contributors-colour](../public_html/assets/img/contributors-colour) directory.
- The `white SVG logo` needs to be added to the [contributors-white](../public_html/assets/img/contributors-white/) directory.

## Contribution workflow

If you'd like to write some code for nf-core/website, the standard workflow
is as follows:

1. Check that there isn't already an issue about your idea in the
   [nf-core/website issues](https://github.com/nf-core/website/issues) to avoid
   duplicating work.
   - If there isn't one already, please create one so that others know you're working on this
2. Fork the [nf-core/website repository](https://github.com/nf-core/website) to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request against the `main` branch and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [basic docs from GitHub](https://help.github.com/articles/fork-a-repo/) or even their [excellent interactive tutorial](https://try.github.io/).

## Tests

When you create a pull request with changes, Github Actions will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.
For now, the only test is for Markdown syntax, using the `markdownlint` package.

## Getting help

For further information/help, please consult the [nf-core/website documentation](https://github.com/nf-core/website#documentation) and don't hesitate to get in touch on the nf-core `tools` channel on [Slack](https://nf-co.re/join/slack/).
