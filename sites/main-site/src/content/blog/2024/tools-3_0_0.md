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

# Avoiding merge conflicts with the new template customisation

If you don't use don't use any of the following template features:

- fastqc
- multiqc
- igenomes
- nf_schema

you can avoid some merge conflicts with a quick update and an intermediate sync:

1. update the template to the latest version.

```bash
nf-core pipelines sync
```

1. Pull the updated `.nf-core.yml` file from the TEMPLATE branch.

```bash
git checkout TEMPLATE -- .nf-core.yml
```

1. add `fastqc`, `igenomes` or `nf_schema` to skip_features.

```yaml title=".nf-core.yml"
template:
  skip_features:
    - fastqc
    - igenomes
    - nf_schema
```

1. Commit the changes

```bash
git add .nf-core.yml
git commit -m "Skip fastqc, igenomes and nf_schema"
```

1. Sync the pipeline again

```bash
run nf-core pipelines sync
```

1. Now you can merge the new template version with _less_ conflicts.

```bash
git merge TEMPLATE
```

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

### How to resolve

Double check that all authors are included.
If that is the case, ignore the incoming changes and keep the credits as they are.

## `.nf-core.yml`

There might be conflicts due to changed order, renamed and new fields, especially in the `template` section.
The following fields were added to the file:

```yaml
bump_version: null
lint: null
nf_core_version: 3.0.0
org_path: null
template:
  author: Author Name
  description: The description of the pipeline
  force: false
  is_nfcore: true
  name: pipelinename
  org: nf-core
  outdir: .
  skip_features: []
  version: 3.0.0
update: null
```

# Relevant template updates

## The `check_max()` function has been removed

The `check_map()` function has been replaced by the core Nextflow funcitonality `resourceLimits`.

The `resourceLimits` are specified in the `nextflow.config` file, and you can remove all mentions to `check_max()` and the parameters that nf-core was using (`max_cpus`, `max_memory` and `max_time`).
You can find more information in the [Nextflow documentation](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits).

## The `nf-validation plugin` has been replaced by `nf-schema`

The `nf-validation` plugin is deprecated and we are now using the new version, the `nf-schema` plugin.

This plugin uses a new JSON schmea draft (2020-12), and thus there are some changes required in the `nextflow_schema.json` and `asses/schema_input.json` files. You can follow the [migration guide](https://nextflow-io.github.io/nf-schema/2.0/migration_guide/) to see the required changes.

As part of these change, the validation parameters are also replaced, the following parameters have been removed from `nextflow.config` and `nextflow_schema.json`:

- validationFailUnrecognisedParams
- validationLenientMode
- validationSchemaIgnoreParams
- validationShowHiddenParams
- validate_params

Instead, the `validation` scope is used to provide `nf-schema` options.

Note that the definition of plugins and the use of `validation` scope has been relocated after the `manifest` scome. This is to allow accessing `manifest` variables to customise the help message.

The use of the new `nf-schema` plugin means that we are also replacing the old `UTILS_NFVALIDATION_PLUGIN` subworkflow by `UTILS_NFSCHEMA_PLUGIN`, and how the input samplesheet is read. Please check [the docs](https://nextflow-io.github.io/nf-schema/2.0/migration_guide/#__tabbed_2_2) to find a description of how to use the new `samplesheetToList()` function.

# Removing for loops and try/catch blocs from `nextflow.config`

To prepare for some soon to come Nextflow changes, code has been reduced from config files.

- Including nf-core configs is not done with a try/catch bloc, and the location these include statements has been moved after defining the profiles. This change is important and it is to make sure that overriding of profiles happens correctly, you can find more information in [this](https://github.com/nextflow-io/nextflow/issues/1792) and [this](https://github.com/nextflow-io/nextflow/issues/5306) Nextflow issue.
