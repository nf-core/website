---
title: "nf-core/tools - 4.0.0"
subtitle: "Friends don't lie (they lint)"
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

This is a major release that requires Nextflows stricter syntax (including a bump of the minimum Nextflow version to 25.10.8), adds more command aliases and removes some old commands (including pytest migration tooling for modules).

As always, if you have any problems or run into any bugs, reach out on the [#tools slack channel](https://nfcore.slack.com/channels/tools).

# Highlights

- [Nextflow strict syntax in the pipeline template](#nextflow-strict-syntax-in-the-pipeline-template)
- [New command aliases](#new-command-aliases)
- Migration to nf-slack

## Webhook notifications removed

The pipeline template no longer includes webhook-based notifications (`hook_url`, `slackreport`, `adaptivecard`).

Instead, use dedicated Nextflow plugins that are purpose-built for this:

- **Slack**: [nf-slack](https://github.com/seqeralabs/nf-slack) (added in the template)
- **Microsoft Teams**: [nf-teams](https://github.com/nvnieuwk/nf-teams)

These plugins integrate directly with Nextflow and give you notification support without coupling it to the pipeline template.

## New command aliases

Sub-subcommands for `modules`, `subworkflows`, and `pipelines` now have shorter aliases so you can type less at the terminal.

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
|                            | `schema lint`                  | `l`                               |
|                            | `sync`                         | `s`                               |
| `test-datasets`            | (command group)                | `t`, `td`, `tds`, `test-datasets` |
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

## Stricter linting for modules

### Linting for `meta` and `ext` keys

The linter now checks that `meta` and `ext` variables are used consistently in module `main.nf` files.
Incorrect or misspelled key names are caught early before they cause hard-to-debug runtime issues.

## Nextflow strict syntax in the pipeline template

Nextflow's [strict syntax mode](https://docs.seqera.io/nextflow/strict-syntax) is becoming the default from Nextflow v26.04.
It enforces more consistent, unambiguous Nextflow code and improves error messages — and nf-core pipelines need to be ready for it.

With this release, the pipeline template is strict syntax compliant and ships with a pre-commit hook that rejects any new code that violates it.
This means issues are caught at commit time rather than failing in CI or, worse, at runtime for your users after they upgrade Nextflow.

If you want to understand what strict syntax means in practice and how to update your pipeline, the [migration guide](/docs/tutorials/migrate_to_strict_syntax/config_strict_syntax) covers the most common patterns.
For the broader picture of how nf-core is rolling out strict syntax across pipelines, modules, and configs, see the [Nextflow syntax roadmap blog post](https://nf-co.re/blog/2025/nextflow_syntax_nf-core_roadmap).

## Other improvements

### Switch from pre-commit to prek in the pipeline template

The pipeline template's GitHub Actions linting workflows now use [prek](https://prek.j178.dev/) instead of `pre-commit` to run git hooks.
`prek` is a fast, Rust-based drop-in replacement for `pre-commit` that uses the same `.pre-commit-config.yaml` format, so no changes are needed to your configuration.
It comes included with nf-core/tools and to replace your pre-commit setup with prek, you run

```bash
prek install --overwrite
```

and you should see moderately to slighlty faster pre-commit checks.

### Sync: respects an existing `defaultBranch`

`nf-core pipelines sync` no longer overwrites the `defaultBranch` setting in `nextflow.config` if it is already present.
This prevents sync from clobbering intentional per-pipeline branch configuration.

# Changelog

You can find the complete changelog and technical details [on GitHub](https://github.com/nf-core/tools/releases/tag/4.0.0).
