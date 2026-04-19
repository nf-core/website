---
title: "Chapter 4: Generating boilerplate files"
subtitle: "How to generate the module file skeleton"
shortTitle: "Chapter 4: Boilerplate files"
---

This chapter explains how to generate an nf-core module template and what each file contains.

## Preparation

Fork and clone the [nf-core/modules GitHub repository](https://github.com/nf-core/modules) to your working environment, then create a new branch for your module:

```bash
git switch master ## ensure you have the latest state of the repository
git switch -c <MY_NEW_MODULE_BRANCH>
```

## Generate the boilerplate files

From the root of the repository, run:

```bash
nf-core modules create <toolname>/<subcommand>
```

### Naming conventions

- All parts of the module name must be lowercase, alphanumeric, with no punctuation or special characters.
- Single-command tools use the tool name only. For a tool executed with `fastp -i <input> -o <output>`, run `nf-core modules create fastp`.
- Tools with subcommands use `<tool>/<subcommand>`, even if you only plan to wrap one subcommand. For `samtools view`, run `nf-core modules create samtools/view`.
- For a third level of subcommand, append it to the subcommand name. For `samtools view flagstats`, run `nf-core modules create samtools/viewflagstat`.

For example, to create a module for the tool `drep` with the subcommand `compare`:

```bash
nf-core modules create drep/compare
```

### Prompts from `nf-core modules create`

The command tries to pre-fill the boilerplate. It searches [Bioconda](https://bioconda.github.io/) and [biocontainers](https://biocontainers.pro/) for the latest version of your tool and adds the container definitions automatically.

You will then be prompted for:

- **Your GitHub username.**
- **A process resource label.** These standardised tags map to default memory, CPU, and wall time. Choose the label that best matches your tool's typical requirements. The defaults for each label are defined in the pipeline template [base.config](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L29-L54). Pipelines can override these values.
- **Whether the module should use a [meta map](https://nf-co.re/docs/developing/components/meta-map).** Answer yes in most cases. Meta maps carry sample metadata alongside files, which pipelines use to drive downstream processing decisions.

:::tip
You can pass any of these values as command-line flags to skip the prompts.
:::

## Output

After the command completes, you will see the following files and directories:

```tree {8-13}
modules/nf-core/drep/compare/
├── environment.yml
├── main.nf
├── meta.yml
└── tests
    └── main.nf.test
```

If you later create a second subcommand (`dereplicate`), the directory structure becomes:

```tree {8-13}
modules/nf-core/drep/
├── compare
│   ├── environment.yml
│   ├── main.nf
│   ├── meta.yml
│   └── tests
│       └── main.nf.test
└── dereplicate
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
```

The next chapter walks through each generated file and explains what to change.
