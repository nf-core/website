---
title: Logo Design Guidelines
subtitle: Guidance on how to use the nf-core logo and variations.
---

All files for the nf-core logo can be seen at the [nf-core/logos GitHub repository](https://github.com/nf-core/logos) and are released under the MIT license.

## nf-core

| Type                         | Logo                                                                                                                                                 | Links                                                                                                                                                                                                                                                                                                                  |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Main nf-core logo            | ![nf-core logo](https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo.png)                                               | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo.ai)                                  |
| Dark backgrounds             | <img class="bg-dark d-block p-2" src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-darkbg.png" width="300">     | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-darkbg.ai)             |
| Monochrome - dark background | <img class="bg-dark d-block p-2" src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-mono-white.png" width="300"> | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-mono-white.ai) |
| Square                       | <img src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-square.png" width="125">                                 | [PNG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.png) <br> [SVG](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.svg) <br> [Adobe Illustrator `.ai`](https://github.com/nf-core/logos/blob/master/nf-core-logos/nf-core-logo-square.ai)             |

## Pipeline logos

### Default

The `nf-core pipelines create` command will generate a logo for your pipeline automatically.

You can manually re-generate different versions of the logo with `nf-core pipelines create-logo`, such as

```bash
nf-core pipelines create-logo --theme dark --format svg isawesome
```

You can also quickly generate a PNG logo yourself by visiting [https://nf-co.re/logo/pipelinename](https://nf-co.re/logo/pipelinename), where `pipelinename` can be anything.

### Custom pipeline logos

Some developers may want to further add branding to their own pipelines by extending the default nf-core logo to include their own mark.

If any any point you want advice, feel free to ask on the nf-core slack under [#graphics](https://nfcore.slack.com/archives/C021QHJRJE8).

If you wish to do so, please use the following guidelines when creating your mark:

1. Simplicity: Keep the mark simple, so easy to visualise
   - A guiding principle in this regard is whether the mark remains recognizable even when presented as an emoji. (128x128px)
2. Colouring: the number of colours should be kept to a minimum, typically following those in the nf-core logo colour scheme
3. Themes: you should generate both light (black text) and dark (white text) mode versions of the logo with mark
4. Positioning: the mark can go in on of the following places
   - to the left of the bottom line of text (pipeline name), with height approximately the same heigh as that line
   - to the left of the entire logo with the height the same of the distance between the top of the top line of text and bottom line of the bottom text (i.e. spanning both the 'nf-core' and pipeline name lines)
5. Other modifications: Do not modify the nf-core apple core in the default logo
6. File formats: you should save your files in the following formats
   - Vector (ideally SVG)
   - Raster (ideally PNG with transparent background at 300 DPI)

If you want to create a 'hex' variant of your logo (e.g. for stickers)

- Make both dark and light version
- Your mark goes above the default nf-core/<pipeline> logo
- The border should ideally be in the nf-core green

Good examples of acceptable pipeline logos:

- [nf-core/sarek](https://raw.githubusercontent.com/nf-core/sarek/master/docs/images/nf-core-sarek_logo_light.png)
- [nf-core/eager](https://raw.githubusercontent.com/nf-core/eager/refs/tags/2.5.2/docs/images/nf-core_eager_logo_outline_drop.png)
- [nf-core/taxprofiler](https://raw.githubusercontent.com/nf-core/taxprofiler/master/docs/images/nf-core-taxprofiler_logo_custom_light.png)
