---
title: "nf-core/tools - 4.0.0"
subtitle: "~Friends~ Linters don't lie"
headerImage: "/assets/images/blog/tools-4_0_0/stricter_syntax.png"
headerImageAlt: "'Stricter Syntax' written in the style of the Stranger Things logo."
pubDate: 2026-04-20T12:00:00+02:00
authors:
  - "mashehu"
label:
  - "tools"
maxHeadingDepth: 3
embedHeaderImage: true
---

This is a major release where we require Nextflow's [strict syntax](https://www.nextflow.io/docs/latest/strict-syntax.html) (including a bump of the minimum Nextflow version to `25.10.4`), prepare for the upcoming [Seqera container transition for modules](/blog/2024/seqera-containers-part-2), add more command aliases, and remove some old commands (including pytest migration tooling for modules).

As always, if you have any problems or run into any bugs, reach out on the [#tools slack channel](https://nfcore.slack.com/channels/tools).

# Highlights

- [Nextflow strict syntax in the pipeline template](#nextflow-strict-syntax-in-the-pipeline-template)
- [Auto-generated container config files](#auto-generated-container-config-files)
- [New command aliases](#new-command-aliases)
- [Webhook notifications removed](#webhook-notifications-removed)

## Nextflow strict syntax in the pipeline template

Nextflow's [strict syntax mode](https://docs.seqera.io/nextflow/strict-syntax) is becoming the default from Nextflow v26.04.
It enforces more consistent, unambiguous Nextflow code and improves error messages, and nf-core pipelines will be ready for it.

With this release, the pipeline template is strict syntax compliant and ships with a pre-commit hook that runs `nextflow lint{:bash}` on your code.

If you want to understand what strict syntax means in practice and how to update your pipeline, the [migration guide](/docs/tutorials/migrate_to_strict_syntax/config_strict_syntax) covers the most common patterns.
For the broader picture of how nf-core is adapting to changes in Nextflow, see the [nf-core syntax adoption roadmap blog post](/blog/2025/nextflow_syntax_nf-core_roadmap).

## Auto-generated container config files

nf-core/tools now generates config files specifying the containers for each platform (Docker, Singularity/Apptainer and Conda) and architecture (AMD64, ARM64) that will be used by Nextflow automatically.

Whenever you install, update, or remove a module, we generate the corresponding container config files based on the container section in the `meta.yml` file (currently only applicable for [FastQC](https://github.com/nf-core/modules/blob/79b36b51048048374b642289bfe9e591ef56fe05/modules/nf-core/fastqc/meta.yml#L88) and [MultiQC](https://github.com/nf-core/modules/blob/008f9d3e61209bf995edac3ba531f54e269e1215/modules/nf-core/multiqc/meta.yml#L110) and will be handled automatically by the new nf-core container commands in the future (see the corresponding [container commands PR](https://github.com/nf-core/tools/pull/3954))).
These files are stored in the `conf/` directory.

```tree
conf/
├── base.config
├── containers_conda_lock_files_amd64.config
├── containers_conda_lock_files_arm64.config
├── containers_docker_amd64.config
├── containers_docker_arm64.config
├── containers_singularity_https_amd64.config
├── containers_singularity_https_arm64.config
├── containers_singularity_oras_amd64.config
├── containers_singularity_oras_arm64.config
...
```

You should not edit these config files manually, and Nextflow will automatically use the correct container for each process based on the platform and architecture.

See our [blog post](/blog/2024/seqera-containers-part-2) for more details about how this works and what will change in the future.

## New command aliases

Continuing from previous releases, where we added `p`, `m`, and `s` aliases for the `pipelines`, `modules`, and `subworkflows` main commands, commands below `modules`, `subworkflows`, and `pipelines` now also have shorter aliases so you need to remember ("was it `remove` or `uninstall`?") and type less.

| Command                    | Subcommand                     | Aliases                           |
| -------------------------- | ------------------------------ | --------------------------------- |
| `modules` & `subworkflows` | `bump-versions` (modules only) | `bump-version`, `bump`, `bv`, `b` |
|                            | `create`                       | `c`                               |
|                            | `install`                      | `add`, `i`                        |
|                            | `lint`                         | `l`                               |
|                            | `list`                         | `ls`                              |
|                            | `patch`                        | `p`                               |
|                            | `remove`                       | `uninstall`, `rm`                 |
|                            | `test`                         | `t`                               |
|                            | `update`                       | `up`, `u`                         |
| `pipelines`                | `bump-version`                 | `bump`, `bv`, `b`                 |
|                            | `create`                       | `c`                               |
|                            | `download`                     | `d`                               |
|                            | `lint`                         | `l`                               |
|                            | `list`                         | `ls`                              |
|                            | `schema lint`                  | `schema l`                        |
|                            | `sync`                         | `s`                               |
| `test-datasets`            | (command group)                | `test-dataset`, `tds`, `td`,`t`   |
|                            | `list`                         | `ls`                              |
|                            | `list-branches`                | `lsb`                             |

So instead of:

```bash
nf-core modules install fastqc
nf-core modules remove samtools/sort
nf-core modules update
```

You can now write:

```bash
nf-core modules add fastqc
nf-core module rm samtools/sort
nf-core m up
```

## Webhook notifications removed

The pipeline template no longer includes webhook-based notifications (`hook_url`, `slackreport`, `adaptivecard`).

Instead, we switched to using dedicated Nextflow plugins for this:

- **Slack**: [nf-slack](https://github.com/seqeralabs/nf-slack) (added in the pipeline template)
- **Microsoft Teams**: [nf-teams](https://github.com/nvnieuwk/nf-teams)

These plugins integrate directly with Nextflow and give you notification support without coupling it to the pipeline template.
Thanks to [@FriederikeHanssen](https://github.com/FriederikeHanssen) for adding this!

## Removed deprecated commands

The `--migrate-pytest{:bash}` flag and all related migration tooling have been removed. This flag was introduced to help transition modules from pytest-based tests to nf-test, and [that migration is now complete](/blog/2025/pytest-to-nf-test-migration).

We now require the `pipelines` prefix for all pipeline commands.
Previously deprecated prefix-free commands, such as `nf-core create{:bash}`, have been removed.

## Other improvements

### Switch from pre-commit to prek in the pipeline template

We switched from [pre-commit](https://pre-commit.com/) to [prek](https://prek.j178.dev/) to run git hooks in the tools repository as well as the pipeline template (and all other nf-core repositories).
prek is a fast, Rust-based drop-in replacement for pre-commit that uses the same `.pre-commit-config.yaml` format, so no changes are needed for your existing configuration.
It comes included with nf-core/tools and to replace your pre-commit setup with prek, you run

```bash
prek install --overwrite
```

and you should see moderately to slightly faster pre-commit checks.

### Stricter module linting for `meta` and `ext` keys and required version topics

The linter now checks that `meta` and `ext` variables are used consistently in module `main.nf` files.
Incorrect or misspelled key names are caught early before they cause hard-to-debug runtime issues.
Thanks to [@mahesh-panchal](https://github.com/mahesh-panchal) for adding this!

We also switched the linting of version topics to be strict, meaning that if [a module or subworkflow does not emit a version topic](https://nf-co.re/docs/developing/migration-guides/update-pipelines), it will fail the linting.

### Apptainer support in module template

The module template now includes Apptainer as a container option alongside Singularity. We [switched already all nf-core modules to this new syntax](https://github.com/nf-core/modules/pull/11260/changes) and new modules created with `nf-core modules create{:bash}` will have Apptainer directives included by default.

### RO-Crate: contributors from `manifest.contributors`

`nf-core pipelines rocrate{:bash}` now reads `manifest.contributors` from `nextflow.config` and maps those entries into the RO-Crate metadata. Pipeline contributors declared in the Nextflow manifest will automatically appear in the generated crate without any extra configuration. Thanks to [@muffato](https://github.com/muffato) for adding this feature, that was on my wishlist for a long time!

### Module linting: compressed file syntax in stubs

The linter now checks that stub blocks use syntax according to the [nf-core module specifications](/docs/specifications/components/modules/general#stub-gzip-files-must-use-echo-and-pipe) for gzipped output files.

### AI and LLM contribution guidelines

Following the [nf-core and AI blog post](/blog/2026/statement-on-ai), we added [guidance to the pipeline template for contributors using AI and LLM tools](https://github.com/nf-core/tools/blob/5b54b0860bfcd92ea34c24489246f087e8cbf132/nf_core/pipeline-template/docs/CONTRIBUTING.md).
It covers how to review, test, and mention AI-assisted contributions.
Thanks [@jfy133](https://github.com/jfy133) for adding this!

### Sync: respects an existing `defaultBranch`

`nf-core pipelines sync{:bash}` no longer overwrites the `defaultBranch` setting in `nextflow.config` if it is already present.
This (hopefully) prevents sync from clobbering intentional per-pipeline branch configuration.

# Changelog

You can find the complete changelog and technical details [on GitHub](https://github.com/nf-core/tools/releases/tag/4.0.0).
