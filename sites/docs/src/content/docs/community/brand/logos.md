---
title: Logos
subtitle: nf-core logos and pipeline logo variations
weight: 4
---

All nf-core logo files are available in the [nf-core/logos GitHub repository](https://github.com/nf-core/logos) under the MIT license.

## nf-core logos

| Type                         | Logo                                                                                                                                                 | Links                                                                                                                                                                                                                                                                                                                  |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Main nf-core logo            | ![nf-core logo](https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo.png)                                               | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.ai)                                  |
| Dark backgrounds             | <img class="bg-dark d-block p-2" src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-darkbg.png" width="300">     | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.ai)             |
| Monochrome - dark background | <img class="bg-dark d-block p-2" src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-mono-white.png" width="300"> | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.ai) |
| Square                       | <img src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-square.png" width="125">                                 | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.ai)             |

## Pipeline logos

Pipeline logos can be created using various tools and formats.

### Default

The `nf-core pipelines create` command automatically generates a pipeline logo.

To manually re-generate different logos, fun the `nf-core pipelines create-logo` command. For example:

```bash
nf-core pipelines create-logo --theme dark --format svg isawesome
```

To quickly generate a PNG logo, visit [https://nf-co.re/logo/pipelinename](https://nf-co.re/logo/pipelinename), where `pipelinename` is the name of your pipeline.

### Custom

Some developers may add branding to their own pipelines by extending the default nf-core logo to include their own mark.

If you wish to do so, please use the following guidelines when creating your mark:

1. Simplicity:
   - Keep the mark simple, so easy to visualise.

     :::tip
     A guiding principle in this regard is whether the mark remains recognizable even when presented as an emoji (128x128px).
     :::

1. Colouring:
   - Limit the number of colours.
   - Use the colour scheme where possible.
1. Themes:
   - Create light (black text) and dark (white text) mode versions of the logo with mark.
1. Positioning:
   - Your mark can be added to:
     - The left of the bottom line of text (pipeline name), with height approximately the same heigh as that line.
     - The left of the entire logo with the height the same of the distance between the top of the top line of text and bottom line of the bottom text (i.e. spanning both the 'nf-core' and pipeline name lines).
1. Other modifications:
   - Do not modify the nf-core apple core in the default logo.
1. File formats:
   - Vector (ideally SVG).
   - Raster (ideally PNG with transparent background at 300 DPI).

:::tip
Ask on the nf-core slack under [#graphics](https://nfcore.slack.com/archives/C021QHJRJE8) for advice.
:::

### Hex

To create a 'hex' variant of your logo (e.g., for stickers):

- Make both dark and light version.
- Your mark goes above the default nf-core/<pipeline> logo.
- The border should ideally be in the nf-core green.

For examples of acceptable pipeline logos, see:

- [nf-core/sarek](https://raw.githubusercontent.com/nf-core/sarek/master/docs/images/nf-core-sarek_logo_light.png)
- [nf-core/eager](https://raw.githubusercontent.com/nf-core/eager/refs/tags/2.5.2/docs/images/nf-core_eager_logo_outline_drop.png)
- [nf-core/taxprofiler](https://raw.githubusercontent.com/nf-core/taxprofiler/master/docs/images/nf-core-taxprofiler_logo_custom_light.png)
