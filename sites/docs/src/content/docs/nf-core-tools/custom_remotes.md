---
title: Remote component repositories
subtitle: Custom remote modules and subworkflows
weight: 50
---

The `modules` and `subworkflows` command groups come with two flags for specifying a custom remote:

- `--git-remote <git remote url>`: Specify the repository from which the modules / subworkflows should be fetched as a git URL. Defaults to the github repository of `nf-core/modules`.
- `--branch <branch name>`: Specify the branch from which the modules / subworkflows should be fetched. Defaults to the default branch of your repository.

For example, if you want to install the `fastqc` module from the repository `nf-core/modules-test` hosted at `gitlab.com`, you can use the following command:

```bash
nf-core modules --git-remote git@gitlab.com:nf-core/modules-test.git install fastqc
```

Or if you want to install the `bam_stats_samtools` subworkflow, you can use the following command:

```bash
nf-core subworkflows --git-remote git@gitlab.com:nf-core/modules-test.git --branch subworkflows install bam_stats_samtools
```

Note that a custom remote must follow a similar directory structure to that of `nf-core/moduleś` for the `nf-core modules` / `nf-core subworkflows` commands to work properly.

The directory where modules / subworkflows are installed will be prompted or obtained from `org_path` in the `.nf-core.yml` file if available. If your modules are located at `modules/my-folder/TOOL/SUBTOOL` your `.nf-core.yml` should have:

```yaml
org_path: my-folder
```

Please avoid installing the same tools from two different remotes, as this can lead to further errors.

The modules commands will during initalisation try to pull changes from the remote repositories. If you want to disable this, for example
due to performance reason or if you want to run the commands offline, you can use the flag `--no-pull`. Note however that the commands will
still need to clone repositories that have previously not been used.

## Private remote repositories

You can use the modules command with private remote repositories. Make sure that your local `git` is correctly configured with your private remote
and then specify the remote the same way you would do with a public remote repository.
