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
