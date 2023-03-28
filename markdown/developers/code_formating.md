---
title: Code formatting
subtitle: Forgot to run Prettier ?
---

# Prettier

You may have stumbled into this recurrent error when attempting a Pull Request:

```
[warn] Code style issues found in 2 files. Forgot to run Prettier?
Error: Process completed with exit code 1.
```

NF-Core community aim for a normalized code style. To do so multiple linting tests are checked and one of them is made with [Prettier](https://prettier.io/).
It checks that the code style is correctly formatted, it supports many languages (thus many different files) and can be integrates with most editors.

# Installation and usage

## Command line

One option is to install prettier in your shell environment.
To do so, you can use the following command with conda:

```bash
conda install prettier
```

and run afterwards the command `prettier -w .`to overwrite the file with a good style or use `prettier -c .` to only check if your files follow the requirements.

## With git as a pre-commit hook

Another solution is to make sure that everytime you commit on a repository prettier is ran on the updated files. To do so, you can add a pre-commit hook that will do that for you.
One way to do that is to add a `.pre-commit-config.yaml` file with the following code:

```yaml
- repo: https://github.com/pre-commit/mirrors-prettier
  rev: '' # Use the sha or tag you want to point at
  hooks:
    - id: prettier
```

## Editor plugins

You can also install prettier as an extension of your favorite [code editor](https://prettier.io/docs/en/editors.html).

### In Visual Studio Code

The extension [Prettier](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode) will do the job for you. It is also available in the pack of useful extension [NF-core](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack).

## During your pull request

If you already did your PR, you can also add a comment with `@nf-core-bot fix linting` in your Pull Request and prettier will be used to apply the required fixes to your code.
