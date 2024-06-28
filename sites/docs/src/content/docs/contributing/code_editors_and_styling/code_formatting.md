---
title: Code formatting
subtitle: Forgot to run Prettier?
---

# Prettier

The nf-core community aims for a normalized code style.
To do so multiple code linting tests are checked and one of them is made with [Prettier](https://prettier.io/).
It checks that the code style is correctly formatted, it supports many languages (thus many different files) and can be integrated with most editors.

You may have stumbled into Prettier for the first time when seeing errors like this, when attempting a Pull Request:

```
[warn] Code style issues found in 2 files. Forgot to run Prettier?
Error: Process completed with exit code 1.
```

:::danger{.fa-brands .fa-youtube title="Bytesize talk"}
Prefer seeing this in action? There is a bytesize talk on this topic that you can watch: [Bytesize 41: Code linting tools](https://nf-co.re/events/2022/bytesize-41-prettier)
:::

## Running with an editor plugin

You can install prettier as an extension of your favorite [code editor](https://prettier.io/docs/en/editors.html) and it will fix formatting issues every time you hit save.

For VSCode there is a [Prettier extension on the VSCode Marketplace](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode). It's also included in the [nf-core VSCode extension pack](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack).

## Running with pre-commit

Another solution is to make sure that every time you commit code on a repository, prettier is run on the updated files.
To do this, you can use [pre-commit](https://pre-commit.com/).
Most nf-core repositories should already come with a `.pre-commit-config.yaml` config file that instructs pre-commit to run prettier:

```yaml
- repo: https://github.com/pre-commit/mirrors-prettier
  rev: '' # Use the sha or tag you want to point at
  hooks:
    - id: prettier
```

To use this, make sure that the pre-commit tool is installed:

```bash
pip install pre-commit
```

and then install the pre-commit hook in the repository working directory:

```bash
pre-commit install
```

After this, every time you run `git commit` the pre-commit tests will run and fix your formatting.

:::tip
If any changes are made, the commit is aborted.
To continue, `git add` the new changes and try again.
:::

## Running manually

If you prefer to be old-school, you can just install prettier in your shell environment.
To do so, you can use the following command with conda:

```bash
conda install prettier
```

Once installed, you can run prettier manually:

```bash
prettier -w .
```

This will overwrite files with formatting fixes
To just check without auto-fixing, use `prettier -c .`

## Running during your pull request

If you already created your pull request, you can add a comment with `@nf-core-bot fix linting` and the kind `@nf-core-bot` will run prettier for you and apply the required fixes to your code.

# Ruff

[Ruff](https://docs.astral.sh/ruff/) is a tool for formatting Python code in a quick and efficient way.
It is designed to be opinionated and enforce consistent standards in whitespace and style.

Several nf-core repositories come with a `pyproject.toml` configuration file that sets the formatting defaults.

### Command line

As with Prettier above, Ruff can be installed on the command line:

```bash
pip install ruff
```

Then run Ruff and it will search recursively:

```bash
ruff format
```

You can provide a file path if you prefer, or use `ruff check` to test without editing files.

### In Visual Studio Code

As with Prettier, there is a [Ruff VSCode plugin](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) that will run Ruff for you whenever you save a file.
