---
title: Install VS Code
subtitle: Learn how to install VS Code and useful extensions
shortTitle: VS Code
weight: 5
---

Visual Studio Code (VS Code) is a popular, free, and open-source code editor.
VS Code has a dedicated official extension that is useful for Nextflow and nf-core pipeline development.

## Installation

VS Code is available for Windows, macOS, and Linux. Download and install it from the [official VS Code website](https://code.visualstudio.com/).

:::note
Windows users should first configure WSL for the optimal development experience. See [Set up a WSL development environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment) for more information.
:::

## Recommended extensions

Extensions enhance VS Code with additional features.

To install an extension:

1. Open VS Code
2. Select the Extensions icon in the left sidebar
   - Alternatively, open the command palette via <kbd>Ctrl/Cmd</kbd>+<kbd>Shift</kbd>+<kbd>P</kbd> and type `Extensions: Install Extension`
3. Search for the extension name
4. Select **Install**

### Nextflow extension

The official [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) provides comprehensive language support for Nextflow development, including:

- Syntax highlighting and validation
- Diagnostics and error detection
- Hover hints and documentation
- Code navigation (go to definition, find references)
- Intelligent code completion
- Code formatting
- Symbol renaming
- Parameter schema validation
- DAG (Directed Acyclic Graph) previews

This extension enforces Nextflow syntax standards and significantly improves your development workflow.

### nf-core extension pack

The [nf-core extension pack](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack) is a community-curated collection of useful tools for pipeline development and documentation:

- **Apptainer/Singularity**: Syntax highlighting for Apptainer/Singularity definition files
- **Docker**: Create, manage, and debug containerized applications
- **EditorConfig**: Maintain consistent coding standards across projects
- **gitignore**: Language support for .gitignore files
- **Markdown Extended**: Enhanced markdown previews with special formatting features
- **Nextflow**: Language support for Nextflow development
- **Prettier**: Code formatter using the Prettier tool
- **Rainbow CSV**: Highlight columns within CSV files with different colors
- **Ruff**: Extremely fast Python linter and code formatter, written in Rust
- **Todo Tree**: Shows TODO, FIXME, and similar comment tags in a tree view
- **YAML**: YAML language support with built-in Kubernetes syntax support

### Remote development pack

For developers working with WSL, SSH, or containers, the [Remote Development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) enables full VS Code functionality in remote environments. It includes:

- **Remote - SSH**: Connect to remote machines via SSH
- **Dev Containers**: Develop inside Docker containers
- **WSL**: Work seamlessly in Windows Subsystem for Linux

## Other editors

While VS Code is widely used, several other editors also support Nextflow development:

### Vim

Vim has official plugin support through the [Nextflow language extension](https://github.com/nextflow-io/nextflow/tree/master/editors/vim). Additional useful packages include:

- **[octo.nvim](https://github.com/pwntester/octo.nvim)**: Edit and review GitHub issues and pull requests directly within the editor
- **[kickstart.nvim](https://github.com/nvim-lua/kickstart.nvim)**: A documented, feature-rich configuration starter for Neovim

### Emacs

An [Emacs mode](https://github.com/Emiller88/nextflow-mode) provides syntax highlighting for Nextflow code. Consider using [Doom Emacs](https://github.com/doomemacs/doomemacs) for a modular editor experience. Complementary packages include:

- **[Magit](https://magit.vc/)**: Git management integrated into Emacs
- **[Forge](https://github.com/magit/forge)**: Extends Magit to work with GitHub and GitLab repositories

### Sublime Text

Sublime Text benefits from a [community-built plugin](https://packagecontrol.io/packages/Nextflow) offering:

- Syntax highlighting compliant with DSL2
- Custom commands with integrations to external services
- Auto-completion and information panels
- Common code snippets

The plugin is available through package control.

## Additional resources

For more information about installing install VS Code and other useful extensions, see:

- Nextflow [environment setup](https://nextflow.io/docs/latest/developer-env.html)
- [Set up Dev Containers in VS Code](./dev-containers.md#set-up-in-visual-studio-code)
