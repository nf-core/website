---
title: 4 - Editing markdown
subtitle: How to develop markdown for the website
---

## How to develop markdown for the website

To develop code for the nf-core website, click the green Gitpod button in upper right of the screen of the [website repository](https://github.com/nf-core/nf-co.re) (or just click the [following link](https://gitpod.io/#https://github.com/nf-core/nf-co.re)).
This will open a Gitpod environment for the website. You'll need to sign in with either Github, GitLab or Bitbucket.

> Some website content is held within other repos (eg. pipeline docs are in pipeline repositories).
> The editing process is the same.

Most of the content on the nf-core website is written in markdown.
You can navigate to the Markdown files by using the explorer panel to the left hand side: ![PNG](/assets/markdown_assets/developers/gitpod/explorer.png)

Then follow `nf-co.re/markdown/developers` to open individual markdown files. Click on a name of the file you want to edit.

Gitpod is able to conveniently render markdown within the editor window.
In the nf-core website repo, open one of the markdown pages.
Next, at the top right of the text editor window, click the preview button:

![PNG](/assets/markdown_assets/developers/gitpod/preview.png)

This should open up in a new window the rendering of the markdown. Now you can edit your raw code and see the changes happening live in the preview. The lines should scroll down at the same time.

## Markdown

Finally, if you are editing Markdown, it is useful to activate an extension for linting.
Linting is an automated process to standardise syntax and formatting of code.
Within nf-core, we use [Prettier](https://prettier.io/) to standardise our markdown.

By adding the Prettier plugin to VSCode, your markdown edits will automatically be run through prettier each time you save.
The plugin may already be installed within your Gitpod environment.
If not, you can add it yourself. Click the extension button on the side panel:

![PNG](/assets/markdown_assets/developers/gitpod/extension.png)

Then choose `Prettier - Code formatter`:
![PNG](/assets/markdown_assets/developers/gitpod/prettier-vscode.png)
