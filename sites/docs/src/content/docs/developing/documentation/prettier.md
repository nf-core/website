---
title: Code formatting
subtitle: Using Prettier and Ruff to maintain consistent code style
shortTitle: Code formatting
---

The nf-core community uses automated code formatting tools to maintain a consistent style across repositories.
Use [Prettier](https://prettier.io/) to check that your documentation is correctly formatted.

Learn how to install and run Prettier to correctly format your documentation.

:::danger{.fa-brands .fa-youtube title="Bytesize talk"}
See [Bytesize 41: Code linting tools](https://nf-co.re/events/2022/bytesize-41-prettier) to see Prettier in action.
:::

## Install Prettier

Install Prettier using conda:

```bash
conda install prettier
```

## Editor plugin

You can install Prettier as an extension for your [code editor](https://prettier.io/docs/en/editors.html).
It will fix formatting issues automatically every time you save a file.

For VS Code, install the [Prettier extension](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode).

:::note
The Prettier extension is included in the [nf-core extension pack](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack).
See [nf-core extension pack](../../get_started/environment-setup/vs-code.md#nf-core-extension-pack) for more information.
:::

## Run Prettier

If you have Prettier installed, you can run it manually to check for formatting issues or apply fixes:

1. Check for changes:

   ```bash
   prettier -c .
   ```

1. Write changes:

   ```bash
   prettier -w .
   ```

## Run Prettier With pre-commit

Configure [pre-commit](https://pre-commit.com/) to run Prettier automatically every time you commit code.
Most nf-core repositories already include a `.pre-commit-config.yaml` file:

```yaml
- repo: https://github.com/pre-commit/mirrors-prettier
  rev: <SHA_OR_TAG> # Use the sha or tag you want to point at
  hooks:
    - id: prettier
```

To run Prettier with pre-commit:

1. Install pre-commit:

   ```bash
   pip install pre-commit
   ```

1. Install the pre-commit hook in your repository:

   ```bash
   pre-commit install
   ```

When you run `git commit`, Prettier runs automatically

:::tip
If changes are made, pre-commit aborts the commit.
To continue, run `git add` with the modified files and try the commit again.
:::

## Run Prettier on a pull request

If you've already created a pull request, add a comment with `@nf-core-bot fix linting`.
The nf-core bot will run Prettier and apply the required fixes.
