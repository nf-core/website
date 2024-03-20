---
title: Basic training to create an nf-core pipeline
subtitle: Adding a nf-core module to your pipeline
---

# Building a pipeline from existing components

Nextflow pipelines can be build in a very modular fashion. In nf-core, we have simple building blocks available: nf-core/modules. Usually, they are wrappers around individual tools. In addition, we have subworkflows: smaller pre-build pipeline chunks. You can think about the modules as Lego bricks and subworkflows as pre-build chunks that can be added to various sets. These components are centrally available for all Nextflow pipelines. To make working with them easy, you can use `nf-core/tools`.

## Identify available nf-core modules

The nf-core pipeline template comes with a few nf-core/modules pre-installed. You can list these with the command below:

```bash
nf-core modules list local
```

```
                                         ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


INFO     Modules installed in '.':

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ Module Name                 ┃ Repository                 ┃ Version SHA                 ┃ Message                    ┃ Date       ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
│ custom/dumpsoftwareversions │ https://github.com/nf-cor… │ 911696ea0b62df80e900ef244d… │ Remove quay from           │ 2023-05-04 │
│                             │                            │                             │ biocontainers (#3380)      │            │
│ fastqc                      │ https://github.com/nf-cor… │ bd8092b67b5103bdd52e300f75… │ Add singularity.registry = │ 2023-07-01 │
│                             │                            │                             │ 'quay.io' for tests        │            │
│                             │                            │                             │ (#3499)                    │            │
│ multiqc                     │ https://github.com/nf-cor… │ 911696ea0b62df80e900ef244d… │ Remove quay from           │ 2023-05-04 │
│                             │                            │                             │ biocontainers (#3380)      │            │
└─────────────────────────────┴────────────────────────────┴─────────────────────────────┴────────────────────────────┴────────────┘

```

These version hashes and repository information for the source of the modules are tracked in the modules.json file in the root of the repo. This file will automatically be updated by nf-core/tools when you create, remove or update modules.

Let’s see if all of our modules are up-to-date:

```bash
nf-core modules update
```

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


? Update all modules or a single named module? All modules
? Do you want to view diffs of the proposed changes? No previews, just update everything
INFO     Updating 'nf-core/custom/dumpsoftwareversions'
INFO     Updating 'nf-core/fastqc'
INFO     Updating 'nf-core/multiqc'
INFO     Updates complete ✨
```

You can list all of the modules available on nf-core/modules via the command below but we have added search functionality to the nf-core website to do this too!

```bash
nf-core modules list remote
```

In addition, all modules are listed on the website: [https://nf-co.re/modules](https://nf-co.re/modules)

## Install a remote nf-core module

To install a remote nf-core module, you can first get information about a tool, including the installation command by executing:

```bash
nf-core modules info salmon/index
```

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re

╭─ Module: salmon/index ───────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ 🌐 Repository: https://github.com/nf-core/modules.git │
│ 🔧 Tools: salmon │
│ 📖 Description: Create index for salmon │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╷ ╷
📥 Inputs │Description │Pattern
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
genome_fasta (file) │Fasta file of the reference genome │
╶────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼───────╴
transcriptome_fasta (file)│Fasta file of the reference transcriptome │
╵ ╵
╷ ╷
📤 Outputs │Description │ Pattern
╺━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
index (directory)│Folder containing the star index files │ salmon
╶───────────────────┼───────────────────────────────────────────────────────────────────────���──────────────────────────┼────────────╴
versions (file) │File containing software versions │versions.yml
╵ ╵

💻 Installation command: nf-core modules install salmon/index

```

:::tip{title="Exercise 4 - Identification of available nf-core modules"}

1. Get information abou the nf-core module `salmon/quant`.
   <details>
      <summary>solution 1</summary>

   ```
   nf-core modules info salmon/quant
   ```

    </details>

2. Is there any version of `salmon/quant` already installed locally?
   <details>
      <summary>solution 2</summary>

   ```
   nf-core modules list local
   ```

   If `salmon/quant` is not listed, there is no local version installed.

      </details>
   :::

The output from the info command will among other things give you the nf-core/tools installation command, lets see what it is doing:

```bash
nf-core modules install salmon/index
```

```

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


╭─ Module: salmon/index  ───────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ 🌐 Repository: https://github.com/nf-core/modules.git                                                                             │
│ 🔧 Tools: salmon                                                                                                                  │
│ 📖 Description: Create index for salmon                                                                                           │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
                             ╷                                                                                              ╷
 📥 Inputs                   │Description                                                                                   │Pattern
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
  genome_fasta  (file)       │Fasta file of the reference genome                                                            │
╶────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼───────╴
  transcriptome_fasta  (file)│Fasta file of the reference transcriptome                                                     │
                             ╵                                                                                              ╵
                    ╷                                                                                                  ╷
 📤 Outputs         │Description                                                                                       │     Pattern
╺━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
  index  (directory)│Folder containing the star index files                                                            │      salmon
╶───────────────────┼───────────────────────────────────────────────────────────────────────���──────────────────────────┼────────────╴
  versions  (file)  │File containing software versions                                                                 │versions.yml
                    ╵                                                                                                  ╵

 💻  Installation command: nf-core modules install salmon/index

gitpod /workspace/basic_training/nf-core-demotest (master) $ nf-core modules install salmon/index

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


INFO     Installing 'salmon/index'
INFO     Use the following statement to include this module:

 include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'
```

The module is now installed into the folder `modules/nf-core`. Now open the file `workflow/demotest.nf`. You will find already several `include` statements there from the installed modules (`MultiQC` and `FastQC`):

```bash title="workflow/demotest.nf"

include { FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
```

Now add the above line underneath it:

```bash title="workflow/demotest.nf"

include { FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'

```

This makes the module now available in the workflow script and it can be called with the right input data.

<!-- TODO/TO DISCUSS here the user now needs to know about how to get their fasta. We could do this here or add a new point for this above -->

We can now call the module in our workflow. Let's place it after FastQC:

```bash title="workflow/demotest.nf"

workflow DEMOTEST {

    ...
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    SALMON_INDEX()
```

Now we are still missing an input for our module. In order to build an index, we require the reference fasta. Luckily, the template pipeline has this already all configured, and we can access it by just using `params.fasta` and `view` it to insppect the channel content. (We will see later how to add more input files.)

```bash
    fasta  = Channel.fromPath(params.fasta)

    fasta.view()

    SALMON_INDEX(
        fasta.map{it -> [id:it.getName(), it]}
    )
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())

```

Now what is happening here:

To pass over our input FastA file, we need to do a small channel manipulation. nf-core/modules typically take the input together with a `meta` map. This is just a hashmap that contains relevant information for the analysis, that should be passed around the pipeline. There are a couple of keys that we share across all modules, such as `id`. So in order, to have a valid input for our module, we just use the fasta file name (`it.getName()`) as our `id`. In addition, we collect the versions of the tools that are run in the module. This will allow us later to track all tools and all versions allow us to generate a report.

How test your pipeline:

```bash
nextflow run main.nf -profile test,docker --outdir results
```

You should now see that `SALMON_INDEX` is run.

(lots of steps missing here)
exercise to add a different module would be nice! => salmon/quant!
comparison to simple nextflow pipeline from the basic Nextflow training would be nice!)

:::tip{title="Exercise 5 - Installing a remote module from nf-core"}

1.  Install the nf-core module `adapterremoval`
    <details>
       <summary>solution 1</summary>

    ```bash
    nf-core modules install adapterremoval
    ```

       </details>

2.  Which file(s) were/are added and what does it / do they do?
    <details>
       <summary>solution 2</summary>

    ```
    Installation added the module directory `/workspace/basic_training/nf-core-demotest/modules/nf-core/adapterremoval`:
    .
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
       ├── main.nf.test
       ├── main.nf.test.snap
       └── tags.yml

    The `test` directory contains all information required to perform basic tests for the module, it rarely needs to be changed. `main.nf` is the main workflow file that contains the module code. All input and output variables of the module are described in the `meta.yml` file, whereas the `environment.yml` file contains the dependancies of the module.
    ```

       </details>

3.  Import the installed `adapterremoval` pipeline into your main workflow.
    <details>
       <summary>solution 3</summary>

    ```bash title="workflows/demotest.nf"
      [...]
       /*
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       */

       include { FASTQC                 } from '../modules/nf-core/fastqc/main'
       include { MULTIQC                } from '../modules/nf-core/multiqc/main'
       include { paramsSummaryMap       } from 'plugin/nf-validation'
       include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
       include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
       include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_demotest_pipeline'
       include { ADAPTERREMOVAL         } from '../modules/nf-core/adapterremoval/main'

      [...]

    ```

       </details>

4.  Call the `ADAPTERREMOVAL` process in your workflow
    <details>
       <summary>solution 4</summary>

    ```bash title="workflows/demotest.nf"
    [...]
    FASTQC (
      ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: ADAPTERREMOVAL
    //
    ADAPTERREMOVAL(

    )
    [...]
    ```

       </details>

5.  Add required parameters for `adapterremoval`to the `ADAPTERREMOVAL` process
    <details>
       <summary>solution 5</summary>

    `adapterremoval` requires three input channels: `meta`, `reads` and `adapterlist`, as outlined in the `meta.yml` of the module. `meta` and `reads` are typically given in one channel as a metamap, whereas the `adapterlist` will be it's own channel for which we should give a path. See here:

    ```bash title="adapterremoval/main.nf"
       [...]
       input:
       tuple val(meta), path(reads)
       path(adapterlist)
       [...]
    ```

    The meta map containing the metadata and the reads can be taken directly from the samplesheet as is the case for FastQC, therefore we can give it the input channel `ch_samplesheet`. The `adapterlist` could either be a fixed path, or a parameter that is given on the command line. For now, we will just add a dummy channel called `adapterlist` assuming that it will be a parameter given in the command line. With this, the new module call for adapterremoval looks as follows:

    ```bash title="workflows/demotest.nf"
    [...]
    //
    // MODULE: ADAPTERREMOVAL
    //
    ADAPTERREMOVAL(
       ch_samplesheet
       params.adapterlist
    )
    [...]
    ```

       </details>

6.  Add the input parameter `adapterlist`
    <details>
       <summary>solution 7</summary>
      In order to use `params.adapterlist` we need to add the parameter to the `nextflow.config`.

    ```bash title="nextflow.config"
       /// Global default params, used in configs
       params {

       /// TODO nf-core: Specify your pipeline's command line flags
       /// Input options
       input                      = null
       adapterlist                = null

      [...]
    ```

    Then use the `nf-core schema build` tool to have the new parameter integrated into `nextflow_schema.json`. The output should look as follows.

    ```
    gitpod /workspace/basic_training/nf-core-demotest (master) $ nf-core schema build

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.13.1 - https://nf-co.re

    INFO [✓] Default parameters match schema validation
    INFO [✓] Pipeline schema looks valid (found 32 params)
    ✨ Found 'params.test' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
    ```

    Select y on the final prompt to launch a web browser to edit your schema graphically.

    </details>

7.  Lint your pipeline
    <details>
       <summary>solution 7</summary>

    ```bash
    nf-core lint
    ```

       </details>

8.  Run the pipeline and inspect the results
    <details>
       <summary>solution 8</summary>

    To run the pipeline, be aware that we now need to specify a file containing the adapters. As such, we create a new file called "adapterlist.txt" and add the adapter sequence "[WE NEED AN ADAPTER SEQUENCE HERE]" to it. Then we can run the pipeline as follows:

    ```bash
    nextflow run nf-core-demotest/ -profile test,docker --outdir test_results --adapterlist /path/to/adapterlist.txt

    ```

       </details>

9.  Commit the changes
    <details>
       <summary>solution 9</summary>

    ```bash
    git add .
    git commit -m "add adapterremoval module"
    ```

       </details>

:::

<p class="text-center">
  <a href="../template_walk_through/" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to Chapter 3
  </a>
  <a href="../add_custom_module/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 5 >
  </a>
</p>
