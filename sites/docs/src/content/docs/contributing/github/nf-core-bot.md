---
title: "🤖 commands for nf-core-bot"
subtitle: "Controlling automated GitHub reactions"
---

# Introduction

To make certain GitHub actions steps a bit easier, there are certain commands that trigger certain code fixes.
You trigger these commands by adding a comment to the Pull Request starting with the specified words.

:::danger
Not all repositories have the same commands available.
:::

# Available commands

## `@nf-core-bot fix linting`

Anytime prettier complains about some incorrect code formatting, you can trigger a fix by adding the following comment to the Pull Request:

```bash
@nf-core-bot fix linting
```

This will start a GitHub Actions workflow that will first react with 👀, run prettier, commit the changes, push them to the repository and then react with 👍🏻 if everything went well.

:::info
Available in [all nf-core pipelines](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/.github/workflows/fix-linting.yml), [nf-core/modules](https://github.com/nf-core/modules/blob/master/.github/workflows/fix-linting.yml), [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/fix-linting.yml) and [nf-core/website](https://github.com/nf-core/website/blob/main/.github/workflows/fix-linting.yml).
:::

## `@nf-core-bot update gpu snapshot path: $PATH`

To update outdated nf-test snapshots for GPU-based tests, you can add the following comment to the Pull Request:

```bash
@nf-core-bot update gpu snapshot path: $PATH
```

where `$PATH` is the path to the test file, e.g. `modules/nf-core/parabricks/applybqsr/tests/main.nf.test`.

This will start a GitHub Actions workflow that will first react with 👀, commit the changes and push them to the repository and then react with 🎉 if everything went well.

:::info
Available in [nf-core/modules](https://github.com/nf-core/modules/blob/master/.github/workflows/update-gpu-snapshot.yml).
:::

## `@nf-core-bot update changelog`

To update the auto-generated changelog entry after you updated the title (or by adding it after the magic words), you can add the following comment to the Pull Request:

```bash
@nf-core-bot update changelog $NEW_TITLE
```

This will update the changelog entry with the title of the Pull Request or the value of `$NEW_TITLE` if provided.

:::info
Available in [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/changelog.yml).
:::

## `@nf-core-bot update textual snapshots`

If the Textual snapshots (run by `tests/pipelines/test_crate_app.py`) fail, an HTML report is generated and uploaded as an artifact.
You can automatically update the snapshots from the PR by posting a comment with the magic words:

```bash
@nf-core-bot update textual snapshots
```

:::warning
Please always check the HTML report to make sure that the changes are expected.
:::

:::info
Available in [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/update-textual-snapshots.yml).
:::
