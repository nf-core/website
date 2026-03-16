---
title: Overview
subtitle: Develop nf-core pipelines
shortTitle: Overview
weight: 1
---

nf-core pipelines use a shared template, a library of reusable components, and community-maintained standards for testing, documentation, and configuration.

This section covers everything you need to start developing best practise pipelines using nf-core.

## Pipelines

nf-core provides a pipeline template and associated tooling to help you build well-structured, community-ready pipelines.
The template enforces consistent file structure, CI/CD workflows, and coding conventions so that all nf-core pipelines behave predictably and are easy to contribute to.

- **[Pipeline file structure](pipelines/template-files):** Understand the files and directories that the nf-core pipeline template generates and their purpose
- **[Adding modules](pipelines/adding-modules):** Integrate standardised nf-core modules into your pipeline using nf-core tools
- **[Release procedure](pipelines/release-procedure):** Step-by-step checklist for preparing and publishing a new pipeline release
- **[External use of nf-core resources](pipelines/external-use):** Guidelines for using nf-core code in non-nf-core pipelines with proper attribution and branding
- **[Renaming branches](pipelines/master_to_main):** Instructions for switching a pipeline's default branch from `master` to `main`

## Components

nf-core modules and subworkflows are reusable building blocks shared across pipelines.
Each module wraps a single tool and follows strict guidelines for metadata, containers, and testing, making it straightforward to combine components into larger workflows.

- **[Creating components](components/creating-components):** Write and contribute reusable modules and subworkflows for use across nf-core pipelines
- **[Storing sample metadata](components/meta-map):** Store and manage sample metadata in nf-core modules using the meta map
- **[ext arguments](components/ext-args):** Use the `ext` directive to inject custom command-line arguments into module scripts
- **[Cross-organisational components](components/cross-org-components):** Combine modules from different organisations within a single subworkflow
- **[Setting custom remotes](components/custom-remotes):** Configure custom git repositories and branches when installing modules and subworkflows
- **[Writing pipeline methods text for MultiQC](components/automated-methods):** Customise automated methods descriptions in MultiQC reports with pipeline and tool citations
- **[Deprecating modules](components/deprecating-components):** Deprecate outdated modules or subworkflows while maintaining backward compatibility

## Containers

nf-core pipelines rely on containers to ensure software reproducibility across different computing environments.
Each module specifies Docker and Singularity containers, which Bioconda and Seqera Containers build and distribute automatically.

- **[Seqera Containers](containers/seqera-containers):** Obtain Docker and Singularity containers from Seqera Containers for use in nf-core modules
- **[ARM64 on Bioconda](containers/arm64-on-bioconda):** Enable ARM64 architecture support in Bioconda and Conda-forge packages for cost-effective cloud instances

## Testing

nf-core uses nf-test to verify that pipelines, modules, and workflows behave as expected.
Tests capture output snapshots and run automatically on every pull request, helping to catch regressions early and giving reviewers confidence that changes are correct.

- **[Overview](testing/overview):** Introduction to the nf-test framework and snapshot-based testing for Nextflow components
- **[Assertions](testing/assertions):** Best practices for writing assertions in nf-test, including snapshot capture and encapsulation requirements
- **[Advanced techniques](testing/advanced):** Advanced assertion patterns for handling complex testing scenarios beyond standard snapshot comparisons

## Template syncs

nf-core pipelines stay aligned with community standards through automated template synchronisation.
When the nf-core template updates, nf-core automatically opens a pull request against each pipeline to propagate improvements, bug fixes, and new CI checks consistently across the community.

- **[Overview](template-syncs/overview):** How the automated template synchronisation system keeps nf-core pipelines up to date
- **[Merging automated PRs](template-syncs/merge-automated-pull-requests):** Review and merge automated template synchronisation pull requests, including conflict resolution
- **[Manually syncing your pipeline](template-syncs/manual-sync):** Trigger manual template synchronisation when automated syncs are unavailable or for custom pipelines
- **[Setting up a pipeline sync retrospectively](template-syncs/set-up-pipeline-sync):** Configure the TEMPLATE branch for pipelines not originally created with `nf-core pipelines create`
- **[Fixing a broken TEMPLATE branch](template-syncs/fix-broken-template-branch):** Resolve issues when the TEMPLATE branch accidentally contains pipeline-specific code

## Documentation

nf-core follows consistent style and formatting standards across all pipelines and components.
Clear, well-structured documentation makes it easier for users to run pipelines correctly and for contributors to understand and extend the code.

- **[Style guide](documentation/style-guide):** Writing standards covering voice, grammar, formatting conventions, and text elements for nf-core documentation
- **[Topic types](documentation/topic-types):** Framework for organising documentation content using the concept, task, reference, and troubleshooting model
- **[Markdown on the nf-core website](documentation/website-markdown):** GitHub Flavored Markdown with additional features that the nf-core website supports
- **[Code formatting](documentation/prettier):** Use Prettier to automatically format documentation and maintain consistent code style
- **[Harshil alignment](documentation/harshil-alignment):** The whitespace-aligned code style used throughout nf-core for consistent Nextflow formatting

## Institutional profiles

Institutional profiles allow sites to share cluster-specific Nextflow configuration with the wider community via nf-core/configs.

- **[Overview](institutional-profiles/overview):** Introduction to institutional profiles as centralised configurations for cluster-specific Nextflow settings
- **[Preparing to write a profile](institutional-profiles/preparing-to-write):** Pre-development checklist including verification of existing profiles and determining configuration scope
- **[File structure requirements](institutional-profiles/file-structure):** Required file locations and structure for institutional profile configuration and documentation
- **[Configuration file components](institutional-profiles/configuration):** Structure cluster-specific Nextflow configuration using the `params` and `process` scopes
- **[Documentation requirements](institutional-profiles/documentation):** Requirements for documenting institutional profiles with cluster information and usage instructions
- **[Testing profiles](institutional-profiles/testing):** Guidelines for thoroughly testing institutional profiles before submitting to nf-core/configs
- **[Using nf-core/configs outside nf-core](institutional-profiles/outside-nf-core):** Integrate nf-core/configs institutional profiles into custom scripts and workflows

## Migration guides

As Nextflow and nf-core evolve, pipelines and components occasionally need updates to adopt new conventions or features.
These guides walk through specific migrations step by step.

- **[Migrating to topic channels](migration-guides/update-pipelines):** Update modules and pipelines to use Nextflow topic channels for version tracking

## External documentation

There is a range of other useful documentation outside of the nf-core initiative, such as for external tooling regularly used by the community.

- **[Nextflow documentation](https://www.nextflow.io/docs/):** Central Nextflow documentation
- **[Nextflow training](https://training.nextflow.io/latest/):** Self-teaching material for learning how to run and write Nextflow code
- **[nf-test documentation](https://www.nf-test.com/):** Documentation for the testing framework used in modules, subworkflows, and pipelines
