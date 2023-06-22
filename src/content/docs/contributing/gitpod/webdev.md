---
title: Editing Markdown
subtitle: How to develop Markdown for the website
weight: 4
parent: gitpod
---

## How to develop Markdown for the website

To develop code for the nf-core website, click the green Gitpod button in upper right of the screen of the [website repository](https://github.com/nf-core/nf-co.re) (or just click the [following link](https://gitpod.io/#https://github.com/nf-core/nf-co.re)).
This will open a Gitpod environment for the website. You'll need to sign in with either Github, GitLab or Bitbucket.

> Some website content is held within other repos (e.g. pipeline docs are in pipeline repositories).
> The editing process is the same.

Most of the content on the nf-core website is written in Markdown.
You can navigate to the Markdown files by using the Explorer panel to the left hand side: ![PNG](@assets/contributing/gitpod/explorer.png)

Then follow `nf-co.re/markdown/developers` to open individual Markdown files. Click on a name of the file you want to edit.

### Previewing the website

When you open Gitpod from the nf-core website (`master`, a branch, or from a pull-request), the environment should be configured to build and serve a preview of the site with your changes.

> This is a great way to preview changes in a pull-request!

When the environment launches, you should see a terminal window at the bottom with the log output from `docker compose` and a _'Simple Browser'_ tab should open in the main view:
![PNG](@assets/contributing/gitpod/website_preview.png)

If you open up a file to edit, you can drag this pane to the left and have both the editor and the browser side by side.
Click refresh in the browser pane to see your edits appear as you type:
![PNG](@assets/contributing/gitpod/website_preview_2.png)

Note that the URL within the _Simple Browser_ can also be opened up directly in a web browser. This can be shared with others and should persist for as long as the Gitpod environment is running.

We use Gitpod prebuilds to make this environment as fast to load as possible.
However, note that you may find some things that do not work in the preview as they do on the main site.
A notable example is that the pipelines, modules and community statistics pages will be empty.
To get this stuff to work, please see the [First-run](https://github.com/nf-core/nf-co.re#first-run) instructions on the website repo readme.

> Remember that Gitpod environments are not secure. Please do not save any secrets in a `config.ini` within Gitpod.

### Previewing markdown

If you prefer, Gitpod is also able to conveniently render Markdown within the Editor window.
Whilst editing some Markdown, at the top right of the text editor window, click the preview button:
![PNG](@assets/contributing/gitpod/preview.png)

This should open up in a new window the rendering of the Markdown code.
Now you can edit your raw code and see the changes happening live in the preview.
The lines should scroll down at the same time.

This method does not require a running server and is typically much faster (and co-scrolls).
However, the preview Markdown rendering is not always identical to the nf-core website.
As such, it's good to use a combination of these two methods when writing new content.

## Markdown

Finally, if you are editing Markdown, it is useful to activate an extension for linting.
Linting is an automated process to standardise syntax and formatting of code.
Within nf-core, we use [Prettier](https://prettier.io/) to standardise our Markdown.

By adding the Prettier plugin to VSCode, your Markdown edits will automatically be run through Prettier each time you save.
The plugin may already be installed within your Gitpod environment.
If not, you can add it yourself. Click the extension button on the side panel:

![PNG](@assets/contributing/gitpod/extension.png)

Then choose `Prettier - Code formatter`:
![PNG](@assets/contributing/gitpod/prettier-vscode.png)
