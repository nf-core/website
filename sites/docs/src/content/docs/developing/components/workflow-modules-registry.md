---
title: Workflow modules and the Nextflow registry
subtitle: Compare Nextflow workflow modules from the Nextflow registry with nf-core subworkflows
shortTitle: Workflow modules and the registry
---

Nextflow 26.04 introduced a module system and a public [Nextflow module registry](https://registry.nextflow.io).
Nextflow 26.10 extended the module system to workflows, so a module can be a standalone process or a standalone workflow.
A module that contains a workflow is called a _workflow module_.

nf-core subworkflows and Nextflow workflow modules are installed, versioned, and laid out on disk differently.
This page describes both and compares them.

:::note
This page describes the module system in Nextflow 26.10.
See the [Nextflow modules documentation](https://docs.seqera.io/nextflow/modules/) for the current behaviour.
:::

## nf-core subworkflows

An nf-core subworkflow is a combination of multiple modules that perform related tasks, stored in the [nf-core/modules](https://github.com/nf-core/modules) repository.
You install a subworkflow with [nf-core/tools](https://github.com/nf-core/tools):

```bash
nf-core subworkflows install fastq_align_star
```

When you install a subworkflow, nf-core/tools:

- Copies the subworkflow, and any subworkflows it uses, into `subworkflows/nf-core/`
- Copies the modules it uses into `modules/nf-core/`
- Copies each component's `tests/` directory along with the component
- Records the Git commit of each installed component in `modules.json`

Each module is installed once per pipeline, so every subworkflow in a pipeline uses the same copy of each module.
You include the subworkflow with a relative path:

```nextflow title="workflows/pipeline.nf"
include { FASTQ_ALIGN_STAR } from '../subworkflows/nf-core/fastq_align_star/main'
```

See [Creating components](creating-components) and the [subworkflow specifications](../../specifications/components/subworkflows/general) for more information.

## Nextflow workflow modules

The Nextflow module registry is a central repository for sharing Nextflow modules, with module discovery, version tracking, and integrity verification.
Nextflow manages registry modules with the [`nextflow module`](https://docs.seqera.io/nextflow/reference/cli/module) command.

A module's `main.nf` defines a single process or a single named workflow.
A module's `meta.yml` file is the module spec.
For a workflow module, the `requires.modules` field of the spec lists the modules that the workflow includes, each pinned to an exact version.

nf-core modules are mirrored to the registry automatically under the `nf-core` namespace.
Most nf-core subworkflows are also available in the registry as workflow modules, without the `subworkflows/` prefix in their name.
For example, the `fastq_align_star` subworkflow is available as `nf-core/fastq_align_star`.
The following excerpt shows fields from the module spec of `nf-core/fastq_align_star` version `0.0.0-6203cd4` in the registry:

```yaml title="meta.yml"
name: "nf-core/fastq_align_star"
version: "0.0.0-6203cd4"
license: "MIT License"
kind: "Workflow"
requires:
  modules:
    - "nf-core/star/align@0.0.0-6c4ed3a"
    - "nf-core/bam_sort_stats_samtools@0.0.0-6203cd4"
```

In the registry version of the subworkflow, `main.nf` includes these dependencies by name, for example `include { STAR_ALIGN } from 'nf-core/star/align'`.

### Install and include a workflow module

To install a workflow module from the registry and use it in your pipeline:

1. Search the registry:

   ```bash
   nextflow module search "star"
   ```

2. View the details of the latest version of a module:

   ```bash
   nextflow module view nf-core/fastq_align_star
   ```

3. Install the module at a specific version:

   ```bash
   nextflow module install nf-core/fastq_align_star -version 0.0.0-6203cd4
   ```

   Nextflow installs the module into the `modules/` directory of your project, and vendors its dependencies under the module's own `modules/` directory.
   It also creates a `.module-info` file that records installation metadata such as the module checksum and registry URL.

4. Include the module by name:

   ```nextflow title="main.nf"
   include { FASTQ_ALIGN_STAR } from 'nf-core/fastq_align_star'
   ```

   Local modules must be included by relative path, for example `./modules/local/fastqc`.

5. Commit the `modules/` directory to your Git repository.

See [Using modules](https://docs.seqera.io/nextflow/modules/using-modules) in the Nextflow documentation for more information.

To list installed modules with their integrity status, run `nextflow module list`.
To remove a module, run `nextflow module remove <namespace>/<name>`.

### Run a workflow module directly

Run a workflow module as a standalone workflow, without a wrapper pipeline:

```bash
nextflow module run myorg/my-workflow --<input_name> <input-value>
```

Each `take:` input of the workflow becomes a parameter of the same name.
Only [typed workflow modules](https://docs.seqera.io/nextflow/workflow-typed) can be run this way.
Typed workflows require the `nextflow.enable.types` feature flag.
As of 15 October 2026, no subworkflow in `nf-core/modules` enables `nextflow.enable.types`.
The registry version of `nf-core/fastq_align_star` doesn't enable it either, so you can't run it with `nextflow module run`.

### Publish a workflow module

To publish your own workflow module to the registry:

1. Create the module:

   ```bash
   nextflow module create myorg/my-workflow -kind Workflow
   ```

   By default, `nextflow module create` scaffolds a module with typed input and output declarations.
   Use `-legacy` to scaffold a module without static typing.

2. Write the workflow in `main.nf`, then generate the module spec from it:

   ```bash
   nextflow module spec myorg/my-workflow
   ```

   For a workflow module, `nextflow module spec` derives `requires.modules` from the `include` declarations in `main.nf`.

3. Validate the module:

   ```bash
   nextflow module validate myorg/my-workflow
   ```

   `nextflow module validate` checks that:
   - The `main.nf`, `meta.yml`, and `README.md` files are present.
   - The module spec contains the required `name`, `version`, `description`, and `license` fields.
   - Every module that `main.nf` includes is a registry module, not a local module.
   - `requires.modules` lists exactly the modules that `main.nf` includes, each pinned to an exact version.

4. Claim a namespace in the Nextflow module registry under **My Resources** > **My Namespaces**.
   The namespace shows as **PENDING REVIEW** until an administrator approves the request.
5. Create an access token under **My Resources** > **My access tokens**.
6. Publish the module:

   ```bash
   export NXF_REGISTRY_TOKEN=<token>
   nextflow module publish myorg/my-workflow
   ```

   You can also set the token with the `registry.apiKey` config option.
   Use `-dry-run` to validate the module structure and metadata without uploading it to the registry.

See [Developing modules](https://docs.seqera.io/nextflow/modules/developing-modules) and [Module registry](https://docs.seqera.io/nextflow/modules/module-registry) in the Nextflow documentation for more information.

## Comparison

The following table compares nf-core subworkflows with Nextflow workflow modules:

| Feature           | nf-core subworkflows                                    | Nextflow workflow modules                                                |
| ----------------- | ------------------------------------------------------- | ------------------------------------------------------------------------ |
| Source            | `nf-core/modules` Git repository, or a custom remote    | Nextflow module registry                                                 |
| Install command   | `nf-core subworkflows install`                          | `nextflow module install`                                                |
| Nextflow version  | Set by the pipeline in `manifest.nextflowVersion`       | 26.10 or later                                                           |
| Install location  | `subworkflows/nf-core/<name>/`                          | `modules/<namespace>/<name>/`                                            |
| Dependencies      | Installed once per pipeline in `modules/nf-core/`       | Vendored under each workflow module's own `modules/` directory           |
| Version tracking  | Git commit SHA of each component in `modules.json`      | `version` in `meta.yml`, plus a checksum in `.module-info`               |
| Version format    | Git commit SHA                                          | Version number, for example `1.0.0`. nf-core mirrors use `0.0.0-<hash>`. |
| Include statement | Relative path                                           | Module name, for example `'nf-core/fastq_align_star'`                    |
| Publishing        | Contributed to `nf-core/modules` through a pull request | Published to a namespace by its members                                  |
| Run standalone    | Not supported by nf-core/tools                          | Yes, with `nextflow module run`, for typed workflow modules              |

Because nf-core/tools installs each module once per pipeline, every subworkflow in a pipeline uses the same version of each module.
Because Nextflow vendors dependencies for each workflow module, two workflow modules can depend on different versions of the same module without conflict.

## Which to use

### nf-core pipeline developers

The nf-core pipeline lint test [`modules_json`](../../nf-core-tools/api_reference/latest/pipeline_lint_tests/modules_json) requires every module installed from `nf-core/modules` to have an entry in `modules.json` with its commit hash.
nf-core/tools creates these entries when it installs a component.
See [Adding modules to pipelines](../pipelines/adding-modules) for more information.

To make a reusable workflow available to all nf-core pipelines, [contribute it as a subworkflow](../../contributing/contribute-components) to `nf-core/modules`.

### Developers of other Nextflow pipelines

You can use either system in your own pipelines:

- The registry lets you include components by name, pins the dependencies of each workflow module to exact versions, and lets you publish workflows under your own namespace. Workflow modules require Nextflow 26.10 or later.
- nf-core/tools tracks installed components in `modules.json`, provides commands to lint, test, and patch components, and can install components from [custom Git remotes](custom-remotes).

See [External use of nf-core resources](../external-use) for attribution guidelines when you use nf-core components outside nf-core.

## Installation and update behaviour

Note the following behaviour when you install modules from the registry:

- **Install a specific version.**
  `nextflow module install` installs the latest version by default.
  Install modules with an explicit `-version` and commit the `modules/` directory to ensure reproducibility.
- **nf-core version numbers don't follow standard semantic versioning.**
  nf-core components in the registry use `0.0.0-<hash>` versions, where the suffix is a short portion of the nf-core component's commit hash.
- **Local changes block updates.**
  Nextflow doesn't update an installed module that has local modifications.
  Use `-force` to overwrite the local changes.
- **Locally modified dependencies block dependency updates.**
  To change the dependency versions of a workflow module, edit `requires.modules` in its `meta.yml` and run `nextflow module install <namespace>/<name> -update-deps`.
  If a vendored dependency has local modifications, Nextflow raises an error instead of overwriting it.
  You can't combine `-update-deps` with `-force`.
- **Both tools install nf-core modules into `modules/nf-core/`.**
  nf-core/tools records the modules it installs in `modules.json`.
  Nextflow records the modules it installs in `.module-info` files.
