---
title: Create a module
subtitle: Create a new module
shortTitle: create
weight: 70
---

This command creates a new nf-core module from the nf-core module template.
This ensures that your module follows the nf-core guidelines.
The template contains extensive `TODO` messages to walk you through the changes you need to make to the template.

You can create a new module using `nf-core modules create`.

This command can be used both when writing a module for the shared [nf-core/modules](https://github.com/nf-core/modules) repository,
and also when creating local modules for a pipeline.

Which type of repository you are working in is detected by the `repository_type` flag in a `.nf-core.yml` file in the root directory,
set to either `pipeline` or `modules`.
The command will automatically look through parent directories for this file to set the root path, so that you can run the command in a subdirectory.
It will start in the current working directory, or whatever is specified with `--dir <directory>`.

The `nf-core modules create` command will prompt you with the relevant questions in order to create all of the necessary module files.

If a bioconda entry is found, the command will populate the containers for you.
If a [bio.tools](https://bio.tools/) entry is found, the command will try to guess the correct input and output channels and [EDAM ontology](https://edamontology.github.io/edam-browser/#topic_0091) entries for you. Make sure to check if the given information is correct.

By following the EDAM ontology URL, you will find a description of that term. Ontology terms can also be browsed in the [EDAM browser](https://edamontology.github.io/edam-browser/#topic_0091).

<!-- RICH-CODEX
working_dir: tmp
timeout: 10
before_command: git clone https://github.com/nf-core/modules.git && cd modules
fake_command: nf-core modules create fastqc --author @nf-core-bot  --label process_low --meta --force
-->

![`cd modules && nf-core modules create fastqc --author @nf-core-bot  --label process_low --meta --force`](../../../../assets/images/tools/nf-core-modules-create.svg)
