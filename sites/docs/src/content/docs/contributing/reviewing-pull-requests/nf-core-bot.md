---
title: "nf-core-bot"
subtitle: "Control automated GitHub reactions"
---

The nf-core-bot can automate code fixes and other actions in pull requests.
You trigger these commands by commenting on the PR with specific keywords.

:::danger
Command availability varies by repository.
Check the repository's `.github/workflows/` directory to see which bot commands are configured.
:::

## Available commands

### `@nf-core-bot fix linting`

When prettier reports code formatting issues, you can automatically fix them by commenting:

```text
@nf-core-bot fix linting
```

This will start a GitHub Actions workflow that will:

1. React with 👀
1. Run prettier
1. Commit changes
1. Push them to the repository
1. React with 👍🏻 if everything goes well

:::info
Available for:

- [All nf-core pipelines](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/.github/workflows/fix_linting.yml)
- [nf-core/modules](https://github.com/nf-core/modules/blob/master/.github/workflows/fix-linting.yml)
- [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/fix-linting.yml)
- [nf-core/website](https://github.com/nf-core/website/blob/main/.github/workflows/fix-linting.yml)
  :::

### `@nf-core-bot update gpu snapshot`

To update outdated nf-test snapshots for GPU-based tests, comment on the PR:

```text
@nf-core-bot update gpu snapshot
```

Without a path, the bot updates the snapshots of the GPU tests related to the changes in the PR.
It finds tests tagged `gpu` or `gpu_highmem`.

To limit the update to specific tests, add a path and, optionally, an extra tag:

```text
@nf-core-bot update gpu snapshot path: $PATH tag: $TAG
```

- `$PATH` is a test file, a directory, or a glob below a component, for example:
  - `modules/nf-core/parabricks/applybqsr/tests/main.nf.test`
  - `modules/nf-core/parabricks/rnafq2bam/tests`
  - `modules/nf-core/parabricks/*`
- `$TAG` is an nf-test tag to run in addition to `gpu`, for example `gpu_highmem`.

The command must be on its own line, in lowercase, with a space after `path:` and `tag:`, and `path:` before `tag:`.
Paths must start with `modules/nf-core/<name>` or `subworkflows/nf-core/<name>`.

The workflow will:

1. React with 👀
1. Run the selected tests with `--update-snapshot` on a GPU runner (a larger one for `gpu_highmem` tests)
1. Commit and push the changed snapshots

If the command is invalid, no tests match, or a step fails, the bot reacts with 👎.
Check the workflow run for details.

:::note
The bot only updates existing snapshot files and only changes snapshots of the tests that were run.
:::

:::info
Available for: [nf-core/modules](https://github.com/nf-core/modules/blob/master/.github/workflows/update-gpu-snapshot.yml).
:::

### `@nf-core-bot update changelog`

To update the auto-generated changelog entry with a new PR title, comment:

```text
@nf-core-bot update changelog $NEW_TITLE
```

This updates the changelog entry using the PR title, or the value of `$NEW_TITLE` if provided.

:::info
Available for [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/changelog.yml).
:::

### `@nf-core-bot update snapshots`

When Textual snapshot tests fail (from `tests/pipelines/test_crate_app.py`), an HTML report is generated and uploaded as an artifact.
To automatically update the snapshots, comment:

```text
@nf-core-bot update snapshots
```

:::warning
Always review the HTML report to verify the changes are expected before merging.
:::

:::info
Available for [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/update-textual-snapshots.yml).
:::
