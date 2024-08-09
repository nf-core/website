---
title: Tool arguments
subtitle: Configure tool containers and arguments
shortTitle: Tool arguments
weight: 3
---

## Modifying tools

Each tool in an nf-core pipeline come preconfigured with a set arguments for an average user.
The arguments are a great place to start and have been tested as a part of the development process.
However, you may want to modify these to fit your own purposes.

**It is unlikely you will need to edit the pipeline code to configure a tool.**

### Docker registries

nf-core pipelines use `quay.io` as the default docker registry for Docker and Podman images.
When specifying a Docker container, it will pull the image from `quay.io` unless a full URI is specified. For example, if the process container is:

- `biocontainers/fastqc:0.11.7--4`

The image will be pulled from quay.io by default, resulting in a full URI of:

- `quay.io/biocontainers/fastqc:0.11.7--4`

If `docker.registry` is specified, it will be used first. For example, if the config value `docker.registry = 'public.ecr.aws'` is specified the image will be pulled from:

- `public.ecr.aws/biocontainers/fastqc:0.11.7--4`

However, the `docker.registry` setting will be ignored if you specify a full URI:

- `docker.io/biocontainers/fastqc:v0.11.9_cv8`

### Tool versions

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of nf-core pipelines uses one container or conda environment per process.

If you need to use a different version of a tool within a pipeline you need to identify the `process` name and override the Nextflow `container` or `conda` definition using the `withName` declaration.

For example, the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline uses a tool called Pangolin that updates an internal database of COVID-19 lineages quite frequently. It doesn't make sense to re-release the `nf-core/viralrecon` pipeline every time a new version of Pangolin is released.

A user can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via a custom configuration file.

1. Check the default version used by the pipeline in the module file for the tool under `modules/nf-core/` directory of the pipeline. For example, for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags) for Docker or [Galaxy Project](https://depot.galaxyproject.org/singularity/) for Singularity
   - Note the container version tag is identical for both container systems, but must include the 'build' ID (e.g.`--pyhdfd78af_1`)
3. Create the custom config accordingly:

   - For Docker:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Singularity:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Conda:

     ```groovy
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.1.17'
         }
     }
     ```



:::warning
Updated tools may come with major changes and may break a pipeline and/or create missing values in MultiQC version tables.
:::

### Tool arguments

You may wish to understand or update which tool arguments a pipeline uses.

You can sometimes find out what parameters are used in a tool in by checking the longer 'help' description of different pipeline parameters, e.g., by pressing the 'help' button next to [this parameter](https://nf-co.re/funcscan/1.0.1/parameters#annotation_bakta_mincontig) in [nf-core/funcscan](https://nf-co.re/funcscan).

There are two main places that a tool can have a tool argument specified:

- The process `script` block
- The `conf/modules.conf` file

Most arguments (both mandatory or optional) are defined in the `conf/modules.conf` file under the `ext.args` entry. Arguments that are defined in the `conf/modules.conf` file can be flexible modified using custom configuration files.

Arguments specified in `ext.args` are then inserted into the module itself via the `$args` variable in the module's bash code

For example, the `-n` parameter could be added to the `BOWTIE_BUILD` process:

```groovy
process {
    withName: BOWTIE_BUILD {
        ext.args = "-n 0.1"
    }
}
```

:::warning
It is recommended to copy and paste existing arguments in a pipelines `conf/modules.config` file to ensure the pipeline can function as expected.
:::
