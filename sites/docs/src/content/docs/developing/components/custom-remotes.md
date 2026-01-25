---
title: Setting custom remotes
subtitle: Set a custom remote for modules and subworkflows
shortTitle: Setting custom remotes
weight: 4
---

The `modules` and `subworkflows` command groups come with two flags for specifying a custom remote:

- `--git-remote <git remote url>`: Specify the repository to fetch modules and subworkflows from as a git URL. Defaults to the GitHub repository of `nf-core/modules`.
- `--branch <branch name>`: Specify the branch to fetch modules and subworkflows from. Defaults to the default branch of your repository.

For example, if you want to install the `fastqc` module from the repository `nf-core/modules-test` hosted at `gitlab.com`, you can use the following command:

```bash
nf-core modules --git-remote git@gitlab.com:nf-core/modules-test.git install fastqc
```

Or if you want to install the `bam_stats_samtools` subworkflow, you can use the following command:

```bash
nf-core subworkflows --git-remote git@gitlab.com:nf-core/modules-test.git --branch subworkflows install bam_stats_samtools
```

A custom remote must follow a similar directory structure to `nf-core/modules` for the `nf-core modules` and `nf-core subworkflows` commands to work properly.

The directory where modules / subworkflows are installed will be prompted or obtained from `org_path` in the `.nf-core.yml` file if available.
If your modules are located at `modules/my-folder/TOOL/SUBTOOL` your `.nf-core.yml` should have:

```yaml
org_path: my-folder
```

Please avoid installing the same tools from two different remotes, as this can lead to further errors.

The modules commands will try to pull changes from remote repositories during initialisation.
If you want to disable this, for example due to performance reasons or to run the commands offline, you can use the `--no-pull` flag.
The commands will still need to clone repositories that have not been previously used.

## Private remote repositories

You can use the modules command with private remote repositories.
Make sure that your local `git` is correctly configured with your private remote and then specify the remote the same way you would with a public remote repository.
