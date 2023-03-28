---
title: Code editor plugins
subtitle: Additions to your coding environment that can help your workflow.
---

# VSCode

[Visual Studio Code](https://code.visualstudio.com/), or VSCode for short, is a code editor redefined and optimized for building and debugging modern web and cloud applications.
VSCode has a huge ecosystem of packages that can be installed to extend functionality.

There is official Nextflow language support with syntax highlighting and auto-completion code snippets.
To use, just search for `Nextflow` in the VSCode Packages tab search bar, or visit
<https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow> (see also the [source code](https://github.com/nextflow-io/vscode-language-nextflow)).

To make it easier to get up and running with the ones that are helpful when working with nf-core pipelines, we have put together an _"extension pack"_ of community favourites.
You can browse, pick and choose the ones you think look good, or can install them all in a single click.

To use, just search for `nf-core-extensionpack` in the VSCode Packages tab search bar, or visit
<https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack>

The extension pack source code can be found on GitHub and is written in super simple syntax - suggestions and improvements welcome!
<https://github.com/nf-core/vscode-extensionpack>

# Sublime Text

[Sublime Text](https://www.sublimetext.com/) is a sophisticated text editor for code, markup and prose. It has a slick user interface, extraordinary features and amazing performance.

Much like VSCode and other code editors, Sublime Text can make use of plugins that extend the native functionality.
[@peterk87](https://github.com/peterk87) has gone the extra mile and built an incredible plugin for Sublime Text users, complete with:

- Nextflow syntax highlighting (DSL2 compliant)
- Custom commands with integration with external services (nf-core, Bioconda, Biocontainers)
- Auto-completion and rich information panels
- Common code snippets

You can get this plugin here: <https://packagecontrol.io/packages/nextflow>

The source code is on GitHub, contributors welcome:
<https://github.com/nf-core/sublime>

# Atom

[Atom](https://atom.io/) is a hackable text editor for the 21st Century, built by GitHub.

The `language-nextflow` Atom package provides syntax highlighting for Nextflow code:
<https://atom.io/packages/language-nextflow>

Other useful packages include:

- `editorconfig`: <https://atom.io/packages/editorconfig>
- `file-icons`: <https://atom.io/packages/file-icons>
- `linter-markdownlint`: <https://atom.io/packages/linter-markdownlint>
- `python-black`: <https://atom.io/packages/python-black>
- `todo-show`: <https://atom.io/packages/todo-show>
- `trailing-spaces`: <https://atom.io/packages/trailing-spaces>

**IN CASE YOU MISSED IT, ATOM HAS BEEN ARCHIVED AND SUPPORT FOR IT WILL BE LIMITED/ABSENT IN THE NEAR FUTURE.**

# Emacs

An Emacs mode written by [@Emiller88](https://github.com/Emiller88) gives Nextflow syntax highlighting:
<https://github.com/Emiller88/nextflow-mode>

If you're looking to get started with Emacs check out [Doom Emacs](https://github.com/hlissner/doom-emacs). If you like modules, it's the editor for you! Check out [DoomCasts: Emacs Doom Screencasts](https://www.youtube.com/playlist?list=PLhXZp00uXBk4np17N39WvB80zgxlZfVwj) for some intros similar to the nf-core bytesize talks.

Other useful packages:

- [Doom Emacs](https://github.com/doomemacs/doomemacs): An Emacs framework for the stubborn martian hacker.
- [`Magit`](https://magit.vc/): A Git Porcelain inside Emacs. [`Forge`](https://magit.vc/manual/forge/) allows you to work with Git forges, such as Github and Gitlab, from the comfort of Magit and the rest of Emacs.

# Vim

[@LukeGoodsell](https://github.com/LukeGoodsell) has put together a Vim plugin that builds on Groovy syntax highlighting to give support for Nextflow `.nf` files: <https://github.com/LukeGoodsell/nextflow-vim>

[@Mxrcon](https://github.com/Mxrcon) has created a fork that supports DSL2 <https://github.com/Mxrcon/nextflow-vim>

Other useful packages:

- [`octo.nvim`](https://github.com/pwntester/octo.nvim): Edit and review GitHub issues and pull requests from the comfort of your favorite editor
- [`kickstart.nvim`](https://github.com/nvim-lua/kickstart.nvim): A small, documented, and featureful neovim starter config
