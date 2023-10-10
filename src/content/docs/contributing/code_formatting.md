---
title: Code formatting
subtitle: Forgot to run Prettier?
---

# Prettier

You may have stumbled into this error when attempting a Pull Request:

```
[warn] Code style issues found in 2 files. Forgot to run Prettier?
Error: Process completed with exit code 1.
```

The nf-core community aims for a normalized code style.
To do so multiple code linting tests are checked and one of them is made with [Prettier](https://prettier.io/).
It checks that the code style is correctly formatted, it supports many languages (thus many different files) and can be integrated with most editors.

> Prefer seeing this in action? There is a bytesize talk on this topic that you can watch: [Bytesize 41: Code linting tools](https://nf-co.re/events/2022/bytesize-41-prettier)

## Installation and usage

### Command line

One option is to install prettier in your shell environment.
To do so, you can use the following command with conda:

```bash
conda install prettier
```

and afterwards run the command `prettier -w .` to overwrite the file with a good style or use `prettier -c .` to only check if your files follow the requirements.

### With git as a pre-commit hook

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

After this, every time you run `git commit ...` the pre-commit tests will run and fix your formatting.

:::note
If any changes are made, the commit is aborted.
To continue, `git add` the new changes and try again.
:::

## Editor plugins

You can also install prettier as an extension of your favorite [code editor](https://prettier.io/docs/en/editors.html).

### In Visual Studio Code

The extension [Prettier](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode) will do the job for you. It is also available in the pack of useful extension [NF-core](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack).

## During your pull request

If you already did your PR, you can also add a comment with `@nf-core-bot fix linting` in your Pull Request and prettier will be used to apply the required fixes to your code.

# Black

[Black](https://black.readthedocs.io/) is a tool for formatting Python code in a quick and efficient way.
It is designed to be opinionated and enforce consistent standards in whitespace and style.

Most nf-core repositories come with a `pyproject.toml` configuration file that sets the formatting defaults.

### Command line

As with Prettier above, Black can be installed on the command line:

```bash
pip install black
```

Then, pass it a path for your files (typically `.` for the current working directory) and it will search recursively:

```bash
black .
```

### In Visual Studio Code

As with Prettier, there is a VSCode plugin that will run Black for you whenever you save a file. See [Black Formatter on the VSCode Marketplace](https://marketplace.visualstudio.com/items?itemName=ms-python.black-formatter) for details.
