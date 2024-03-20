---
title: Basic training to create an nf-core pipeline
subtitle: Setting up the gitpod environment for the course
---

### Gitpod

For this tutorial we will use Gitpod, which runs in the learners web browser. The Gitpod environment contains a preconfigured Nextflow development environment
which includes a terminal, file editor, file browser, Nextflow, and nf-core tools. To use Gitpod, you will need:

- A GitHub account
- Web browser (Google Chrome, Firefox)
- Internet connection

Click the link and log in using your GitHub account to start the tutorial:

<p class="text-center">
  <a href="https://www.gitpod.io/#https://github.com/nf-core/basic_training" class="btn btn-lg btn-success" target="_blank">
    Launch GitPod
  </a>
</p>

For more information about Gitpod, including how to make your own Gitpod environement, see the Gitpod bytesize talk on youtube (link to the bytesize talk),
check the [nf-core Gitpod documentation](gitpod/index) or [Gitpod's own documentation](https://www.gitpod.io/docs).

<details>
<summary> Expand this section for instructions to explore your Gitpod environment</summary>

#### Explore your Gitpod interface

You should now see something similar to the following:

(insert Gitpod welcome image)

- **The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). Click the Explorer button to see which files are in this repository.
- **The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.
- **The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. You should also see the nf-training material browser (<https://training.nextflow.io/>).

To test that the environment is working correctly, type the following into the terminal:

```bash
nextflow info
```

This should come up with the Nextflow version and runtime information:

```
Version: 23.10.0 build 5889
Created: 15-10-2023 15:07 UTC (15:07 GMT)
System: Linux 6.1.54-060154-generic
Runtime: Groovy 3.0.19 on OpenJDK 64-Bit Server VM 17.0.8-internal+0-adhoc..src
Encoding: UTF-8 (UTF-8)
```

#### Reopening a Gitpod session

When a Gitpod session is not used for a while, i.e., goes idle, it will timeout and close the interface.
You can reopen the environment from <https://gitpod.io/workspaces>. Find your previous environment in the list, then select the ellipsis (three dots icon) and select Open.

If you have saved the URL for your previous Gitpod environment, you can simply open it in your browser.

Alternatively, you can start a new workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

If you have lost your environment, you can find the main scripts used in this tutorial in the `nf-training` directory.

#### Saving files from Gitpod to your local machine

To save any file locally from the explorer panel, right-click the file and select Download.

</details>
<p class="text-center">
  <a href="../" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to training index
  </a>
  <a href="../nf_core_create_tool/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 2 >
  </a>
</p>
