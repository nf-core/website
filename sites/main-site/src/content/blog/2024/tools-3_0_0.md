---
title: 'nf-core/tools - 3.0.0'
subtitle: 'What flavour of the nf-core template do you prefer?'
headerImage: 'https://images.unsplash.com/photo-1583258292688-d0213dc5a3a8'
headerImageAlt: 'Yellow and red apples on black plastic crate'
pubDate: 2024-10-03T12:00:00+01:00
authors:
  - 'mirpedrol'
  - 'mashehu'
label:
  - 'tools'
---

# Highlights

## ⛓️‍💥 Breaking changes

- All pipeline commands now need the `pipelines` prefix. For example, `nf-core lint{:bash}` is now `nf-core pipelines lint{:bash}`. This makes the commands more consistent with `nf-core modules{:bash}` and `nf-core subworkflows{:bash}` commands.

## ✨ New features

- More customisation for pipeline templates. The template has been divided into features which can be skipped, e.g. you can create a new pipeline without any traces of FastQC in it. You can now strip down the pipeline to the bare minimum and add only the tools you need. For nf-core pipelines we still require some core features (e.g. documentation, CI tests, etc.) to be present, but you can still customise the pipeline to your needs.
- A new Text User Interface app when running nf-core pipelines create to help us guide you through the process better (no worries, you can still use the CLI if you give all values as parameters).

- nf-validation replaced nf-schema in the pipeline template.
- CI tests now lint with the nf-core tools version matching the template version of the pipeline, to minimise errors in opened PRs with every new tools release.

## 🫡 Deprecations

- The `nf-core licences{:bash}` command is deprecated.

# Common merge conflicts and how to resolve them

## `README.md`

In the `## Credits` section, you might encounter a merge conflict like this:

```diff
<<<<<<< HEAD
- nf-core/mag was written by [Hadrien Gourlé](https://hadriengourle.com) at [SLU](https://slu.se), [Daniel Straub](https://github.com/d4straub) and - [Sabrina Krakau](https://github.com/skrakau) at the [Quantitative Biology Center (QBiC)](http://qbic.life). [James A. Fellows Yates](https://github.- com/jfy133) and [Maxime Borry](https://github.com/maxibor) at the [Max Planck Institute for Evolutionary Anthropology](https://www.eva.mpg.de) joined in version 2.2.0.
-
- Other code contributors include:
-
- - [Antonia Schuster](https://github.com/AntoniaSchuster)
- - [Alexander Ramos](https://github.com/alxndrdiaz)
- - [Carson Miller](https://github.com/CarsonJM)
- - [Daniel Lundin](https://github.com/erikrikarddaniel)
- - [Danielle Callan](https://github.com/d-callan)
- - [Gregory Sprenger](https://github.com/gregorysprenger)
- - [Jim Downie](https://github.com/prototaxites)
- - [Phil Palmer](https://github.com/PhilPalmer)
- - [@willros](https://github.com/willros)
-
- Long read processing was inspired by [caspargross/HybridAssembly](https://github.com/caspargross/HybridAssembly) written by Caspar Gross [@caspargross](https://github.com/caspargross)
=======
+ nf-core/mag was originally written by Hadrien Gourlé, Daniel Straub, Sabrina Krakau, James A. Fellows Yates, Maxime Borry.
>>>>>>> TEMPLATE
```

just ignore the incoming changes and keep the credits as they are.

## `.nf-core.yml`

There might be conflicts due to changed order, renamed and new fields, especially in the `template` section. In general, you can accept all incoming changes.
The following fields were added to the file:

```yaml
template:
  force: false
  is_nfcore: true
  org: nf-core
  outdir: .
  skip_features: []
  version: 3.0.0
```
