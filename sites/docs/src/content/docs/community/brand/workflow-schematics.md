---
title: Workflow schematics
subtitle: Components for making pipeline diagrams
weight: 5
---

Simple workflow schematics help outline the main functionality and steps of a pipeline.

Most workflow schematics are made with vector image editors, such as the open-source tool [Inkscape](https://inkscape.org/) or commercial suite [Adobe Illustrator](https://www.adobe.com/products/illustrator.html). Useful tools for collaborative prototyping include [Google Drawings](https://docs.google.com/drawings/) and [LucidChart](https://www.lucidchart.com/pages/).

All examples and components provided below can be opened in these editors. Various parts can be borrowed and/or modified. Components are also available on [bioicons](https://bioicons.com/icons/cc-0/Chemo-_and_Bioinformatics/James-A--Fellows-Yates/metromap_style_pipeline_workflow_components.svg), which have direct import extensions for [Inkscape](https://inkscape.org/) and [draw.io](https://app.diagrams.net/).

## Examples

See below for examples of nf-core workflow schematics that can be re-used and modified for your own pipeline.

:::warning
Check for any attributions to be included within any derivative images, as defined by the corresponding license.
:::

Select the schematic image to see the original.

|                                                             Workflow Example                                                             | nf-core Pipeline                                  | License/Publication                                                                                                                                                      | Suggested attribution                           |
| :--------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------- |
|                 ![nf-core/sarek](https://raw.githubusercontent.com/nf-core/sarek/master/docs/images/sarek_workflow.png)                  | [nf-core/sarek](https://nf-co.re/sarek)           | From [Garcia _et al._ (2020, F1000 Research)](https://doi.org/10.12688/f1000research.16665.1) under a [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. |                                                 |
|      ![nf-core/eager workflow simple](https://raw.githubusercontent.com/nf-core/eager/master/docs/images/usage/eager2_workflow.png)      | [nf-core/eager](https://nf-co.re/eager)           | From [Fellows Yates _et al._ (2021, PeerJ)](https://doi.org/10.7717/peerj.10947) under a [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license               | CC-BY 4.0. Design originally by Zandra Fagernäs |
| ![nf-core/eager workflow detailed](https://raw.githubusercontent.com/nf-core/eager/master/docs/images/usage/eager2_metromap_complex.png) | [nf-core/eager](https://nf-co.re/eager)           | From [Fellows Yates _et al._ (2021, PeerJ)](https://doi.org/10.7717/peerj.10947) under a [MIT](https://github.com/nf-core/dualrnaseq/blob/master/LICENSE) license        |                                                 |
| ![nf-core/dualrnaseq workflow](https://raw.githubusercontent.com/nf-core/dualrnaseq/master/docs/images/Workflow_diagram_dualrnaseq.png)  | [nf-core/dualrnaseq](https://nf-co.re/dualrnaseq) | By Regan Hayward under [MIT](https://github.com/nf-core/dualrnaseq/blob/master/LICENSE) license                                                                          |                                                 |
|            ![nf-core/circrna workflow](https://raw.githubusercontent.com/nf-core/circrna/master_old/docs/images/workflow.png)            | [nf-core/circrna](https://nf-co.re/circrna)       | By Barry Digby under [MIT](https://github.com/nf-core/circrna/blob/master/LICENSE) license                                                                               |                                                 |
|             !["nf-core/mag workflow"](https://raw.githubusercontent.com/nf-core/mag/main/docs/images/mag_metromap_dark.png)              | [nf-core/mag](https://nf-co.re/mag)               | By Sabrina Krakau under [MIT](https://github.com/nf-core/mag/blob/master/LICENSE) license                                                                                | CC-BY 4.0. Design originally by Zandra Fagernäs |
|                 ![nf-core/bactmap workflow](https://github.com/nf-core/bactmap/raw/dev/docs/images/bactmap_pipeline.png)                 | [nf-core/bactmap](https://nf-co.re/bactmap)       | By Anthony Underwood under [MIT](https://github.com/nf-core/mag/blob/master/LICENSE) license                                                                             |                                                 |
|    ![nf-core/cutandrun workflow](https://raw.githubusercontent.com/nf-core/cutandrun/3.1/docs/images/cutandrun-flow-diagram-v3.0.png)    | [nf-core/cutandrun](https://nf-co.re/cutandrun)   | By Chris Cheshire under [MIT](https://github.com/nf-core/cutandrun/blob/master/LICENSE) license                                                                          |                                                 |
|             ![nf-core/cutandrun workflow](https://raw.githubusercontent.com/nf-core/sarek/dev/docs/images/sarek_subway.png)              | [nf-core/sarek](https://nf-co.re/sarek)           | By Maxime U Garcia under [MIT](https://github.com/nf-core/sarek/blob/master/LICENSE) license                                                                             |                                                 |
| ![nf-core/rnaseq metro map grey](https://raw.githubusercontent.com/nf-core/rnaseq/master/docs/images/nf-core-rnaseq_metro_map_grey.png)  | [nf-core/rnaseq](https://nf-co.re/rnaseq)         | By Sarah Guinchard under [MIT](https://github.com/nf-core/sarek/blob/master/LICENSE) license                                                                             |                                                 |
|        ![nf-core/isoseq metro map](https://raw.githubusercontent.com/nf-core/isoseq/1.1.4/docs/images/Isoseq_pipeline_metro.png)         | [nf-core/isoseq](https://nf-co.re/isoseq)         | By Sébastien Guizard under [MIT](https://github.com/nf-core/isoseq/blob/master/LICENSE) license                                                                          |                                                 |

## Components

<!-- TODO: Add components table -->

## Use draw.io

The web app [draw.io](https://app.diagrams.net/) helps you create, render and export different diagrams including metro-maps.
For even more convenience, you can use the asset library [nf-core xml item library](https://raw.githubusercontent.com/nf-core/website/refs/heads/main/sites/docs/src/assets/images/graphic_design_assets/workflow_schematics_components/generic/nf-core_components.xml). It contains all of the components from the components above.
To import it to draw.io, select **File > Open library from > URL* and paste:

```bash
https://raw.githubusercontent.com/nf-core/website/refs/heads/main/sites/docs/src/assets/images/graphic_design_assets/workflow_schematics_components/generic/nf-core_components.xml
```

:::tip
Components can also be accessed via [bioicons](https://bioicons.com/icons/cc-0/Chemo-_and_Bioinformatics/James-A--Fellows-Yates/metromap_style_pipeline_workflow_components.svg).
:::
