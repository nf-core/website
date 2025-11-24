---
title: "nf-core-bot"
subtitle: "Learn how to control automated GitHub reactions"
weight: 1
---

The nf-core-bot can automate code fixes and other actions in pull requests.
You trigger these commands by commenting on the PR with specific keywords.

:::danger
Command availability varies by repository. Check the repository's `.github/workflows/` directory to see which bot commands are configured.
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
- [All nf-core pipelines](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/.github/workflows/fix-linting.yml)
- [nf-core/modules](https://github.com/nf-core/modules/blob/master/.github/workflows/fix-linting.yml)
- [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/fix-linting.yml)
- [nf-core/website](https://github.com/nf-core/website/blob/main/.github/workflows/fix-linting.yml)
:::

### `@nf-core-bot update gpu snapshot path: $PATH`

To update outdated nf-test snapshots for GPU-based tests, comment on the PR:

```text
@nf-core-bot update gpu snapshot path: $PATH
```

Replace `$PATH` with the path to the test file. For example, `modules/nf-core/parabricks/applybqsr/tests/main.nf.test`.

The workflow will:

1. React with 👀
1. Update the snapshots
1. Commit and push the changes
1. React with 🎉 on success

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

When Textual snapshot tests fail (from `tests/pipelines/test_crate_app.py`), an HTML report is generated and uploaded as an artifact. To automatically update the snapshots, comment:

```text
@nf-core-bot update snapshots
```

:::warning
Always review the HTML report to verify the changes are expected before merging.
:::

:::info
Available for [nf-core/tools](https://github.com/nf-core/tools/blob/main/.github/workflows/update-textual-snapshots.yml).
:::
