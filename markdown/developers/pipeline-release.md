---
title: Pipeline Release Instructions
subtitle: Instructions for releasing a sanger-tol pipeline
---

## nf-core Release Checklist

This page is heavily inspired by the [nf-core release checklist](https://nf-co.re/docs/contributing/release_checklist).

## Branch model

First a reminder about how we use branches in sanger-tol.

`main` is the release branch. Releases can only happen from the `main` branch
and code present on the `main` branch _has_ to be released. `main` can be
updated by merging either the staging branch, `dev`, or a bugfix branch.

`dev` is the staging branch, which accumulates new features before release.
Code on the `dev` branch should always pass tests and be functional (this is a
condition of our [review guidelines](/docs/contributing/review_checklist)).

Bugfix branches can be used to patch a released pipeline. The branch needs to
be created off `main` and merged into `main` for immediate release. There needs
to be extra caution when merging a bugfix branch into `main` as there isn't the
`dev` branch as a buffer.

Our model only supports 1 "active" release at a time. For instance, if you have
already released 1.0.0 and 1.1.0, 1.1.0 is considered the "active" version, and
the next version, according to [Semantic Versioning](https://semver.org/), can
only be 1.1.1, 1.2.0, or 2.0.0. In this example, releasing 1.0.1 is **not**
supported.

## Release steps

Follow the "Before you release" and "Steps to release" from the [nf-core release checklist](https://nf-co.re/docs/contributing/release_checklist) with the following adaptations:

1. Replace "nf-core" with "sanger-tol".
2. Since we develop directly off the sanger-tol repository, the version bump happens there, not "on your fork".
3. The "core team member" who can activate the integration with Zenodo is [@muffato](https://github.com/muffato).
4. We don't have Twitter integration.

About version numbers and release names:

1. The release _number_ is made of digits and dots **only**, and follows [Semantic Versioning](https://semver.org/).
   It is used as the git _tag_ too.
2. The release _name_ is whatever you want. It can be Harry Potter, Stargate, or Pokemon. Express your creativity! You only need a new name for major and minor releases. Patch releases reuse the last name with the addition of `(patch ${PATCH_NUMBER})`.
3. The release _title_ you enter on GitHub must be of the form `v${RELEASE_NUMBER} - ${RELEASE_NAME}`. No need to include the pipeline name there, especially as it is added by Zenodo later.
4. The release _description_ you enter on GitHub should be the same from what you have in the Changelog, e.g.:

```text
## [[${RELEASE_NUMBER}](https://github.com/sanger-tol/${PIPELINE_NAME}/releases/tag/${RELEASE_NUMBER})] - ${RELEASE_NAME} - [${RELEASE_DATE}]

*Summary of the release*

### Enhancements & fixes

- *List of what's changed. Indicate whether they're bug fixes, additions, or breaking changes*

### Parameters

| Old parameter | New parameter |
| ------------- | ------------- |
|               | --added       |
| --removed     |               |
| --old         | --new         |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency   | Old version | New version |
| ------------ | ----------- | ----------- |
| name_removed | 2.30.0      |             |
| name_added   |             | 5.4.3       |
| name_changed | 0.8.10      | 0.8.11      |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.
```

Alternatively, GitHub can also generate release notes from the list of commits. Just make sure you have this at the top:

```text
## [[${RELEASE_NUMBER}](https://github.com/sanger-tol/${PIPELINE_NAME}/releases/tag/${RELEASE_NUMBER})] - ${RELEASE_NAME} - [${RELEASE_DATE}]

*Summary of the release*
```

## Integration with Zenodo

Once you've made a release, a record is automatically created on Zenodo.
Tell [@muffato](https://github.com/muffato) or [@DLBPointon](https://github.com/DLBPointon) who can then do the following:

1. Check the release notes. Sometimes the conversion from Markdown doesn't work well, especially the tables. If this happens, copy-paste the rendered Markdwon from GitHub into the Zenodo editor.
2. Change the record type from "Software" to "Workflow".
3. Check that all authors are properly named (first name and last name identified), have an ORCiD, and have the recognised Sanger affiliation (it is named "Wellcome Sanger Institute (WTSI)", not "Wellcome Sanger Institute", and has the [ROR](https://ror.org/05cy4wa09) logo).
4. Add the pipeline to the [sanger-tol Zenodo community](https://zenodo.org/communities/sanger-tol).
5. Check that the licence is correctly set to MIT.
6. In the "Software" section, link to GitHub URL, enter "Nextflow" as the language, and set the "Development Status" to "Active".

## Integration with WorkflowHub

As part of an EBP and ERGA recommendation, we should deposit our workflows into [WorkflowHub](https://workflowhub.eu/programmes/37) too.
There isn't an automated way of doing that yet, so in the meantime we need to manually upload the releases.

In order to upload a record, you will need to create an account and ask to join one of our two teams (["Genome Assembly"](https://workflowhub.eu/projects/204) and ["Genome Analysis"](https://workflowhub.eu/projects/205)).
Instructions are on their website: <https://about.workflowhub.eu/docs/>
