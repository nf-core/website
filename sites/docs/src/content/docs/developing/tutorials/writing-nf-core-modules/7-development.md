---
title: "Chapter 7: Development workflow"
subtitle: "How to contribute to the community"
shortTitle: "Chapter 7: Development workflow"
---

This chapter outlines the workflow for contributing to the [nf-core/modules repository](https://github.com/nf-core/modules/).

This training used a fork of that repository, but you can apply everything you have learnt to your own private module repository or to local modules inside a pipeline. The nf-core community encourages you to contribute modules upstream so others can benefit.

## Check

Before writing a module, check that someone has not already built it or started on it:

1. Search the [modules page](https://nf-co.re/modules) on the nf-core website.
2. Search [open pull requests](https://github.com/nf-core/modules/pulls) on the nf-core/modules repository.
3. Search [open issues](https://github.com/nf-core/modules/issues) on the nf-core/modules repository.

Based on what you find:

- **Module exists**: skip the writing and install it into your pipeline.
- **Open PR**: offer to help out, or take over if development has stalled.
- **Open issue**: assign yourself, or ask the current assignee if they want help.
- **Nothing yet**: [open an issue](https://github.com/nf-core/modules/issues) to signal you are working on it.

Once you claim a module, fork and clone the repository, create a new branch, and follow the previous chapters to write and test the module.

## Write

Generate the boilerplate template files:

```bash
nf-core modules create <toolname>/<subcommand>
```

Fill in each file by following the `TODO` comments, the [nf-core module specifications](https://nf-co.re/docs/specifications/components/modules/general), and this training material.

## Test

Run the tests and generate the snapshot:

```bash
nf-core modules test <toolname>/<subcommand>
```

Check the snapshot file looks correct before moving on.

## Lint

Standardisation is central to nf-core. The linting tool flags issues before submission:

```bash
nf-core modules lint <toolname>/<subcommand>
```

The console output lists any issues to fix. Some issues can be auto-fixed with `--fix`, in which case the linter prints the command to run.

Check each fix carefully. The linter may add empty sections that you still need to populate.

:::note
Linting is optional for private custom modules, but still recommended to catch small issues.
:::

:::warning
The linter cannot check every aspect of the [nf-core module specification](https://nf-co.re/docs/specifications/components/modules/general). Review the specification yourself before submission.
:::

## Submit

Commit and push your module to your fork:

```bash
git commit -am "Add new module <toolname>/<subcommand>"
git push
```

Open a [pull request](https://github.com/nf-core/modules/pulls) from your fork and branch to the nf-core/modules repository.

To attract a reviewer, ask in the nf-core Slack [#request-review](https://nfcore.slack.com/archives/CQY2U5QU9) channel. See the [join instructions](https://nf-co.re/join) if you are not already on Slack.

If a reviewer asks for changes, do not close the PR. Push updates to the same branch, then request another review.

Once the PR is approved, merge it and the module is available to everyone.

If you get stuck, ask on the [#modules](https://nfcore.slack.com/archives/CJRH30T6V) channel or the [#nostupidquestions](https://nfcore.slack.com/archives/C043FMKUNLB) channel.

Thank you for contributing to the community.

The final chapter covers how to use modules in your own pipelines.
