---
title: Basic training to create an nf-core pipeline
subtitle: Adding a local custom module to your pipeline
---

### Adding a local module

If there is no nf-core module available for the software you want to include, the nf-core tools package can also aid in the generation of a local module that is specific for your pipeline. To add a local module run the following:

```

nf-core modules create

```

Open ./modules/local/demo/module.nf and start customising this to your needs whilst working your way through the extensive TODO comments! For further help and guidelines for the modules code, check out the [modules specific documentation](https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial).

### Making a local module for a custom script

To generate a module for a custom script you need to follow the same steps when adding a remote module.
Then, you can supply the command for your script in the `script` block but your script needs to be present
and _executable_ in the `bin`
folder of the pipeline.
In the nf-core pipelines,
this folder is in the main directory and you can see in [`rnaseq`](https://github.com/nf-core/rnaseq).
Let's look at an publicly available example in this pipeline,
for instance [`tximport.r`](https://github.com/nf-core/rnaseq/blob/master/bin/tximport.r).
This is an Rscript present in the [`bin`](https://github.com/nf-core/rnaseq/tree/master/bin) of the pipeline.
We can find the module that runs this script in
[`modules/local/tximport`](https://github.com/nf-core/rnaseq/blob/master/modules/local/tximport/main.nf).
As we can see the script is being called in the `script` block, note that `tximport.r` is
being executed as if it was called from the command line and therefore needs to be _executable_.

<blockquote style="border-left: 4px solid #F0AD4E; background-color: #FFF3CD; padding: 10px;">

<h4 style="margin-top: 0;">TL;DR</h4>

1. Write your script on any language (python, bash, R,
   ruby). E.g. `maf2bed.py`
2. If not there yet, move your script to `bin` folder of
   the pipeline and make it
   executable (`chmod +x <filename>`)
3. Create a module with a single process to call your script from within the workflow. E.g. `./modules/local/convert_maf2bed/main.nf`
4. Include your new module in your workflow with the command `include {CONVERT_MAF2BED} from './modules/local/convert_maf2bed/main'` that is written before the workflow call.
</blockquote>

_Tip: Try to follow best practices when writing a script for
reproducibility and maintenance purposes: add the
shebang (e.g. `#!/usr/bin/env python`), and a header
with description and type of license._

### 1. Write your script

Let's create a simple custom script that converts a MAF file to a BED file called `maf2bed.py` and place it in the bin directory of our nf-core-testpipeline::

```

#!/usr/bin/env python
"""bash title="maf2bed.py"
Author: Raquel Manzano - @RaqManzano
Script: Convert MAF to BED format keeping ref and alt info
License: MIT
"""
import argparse
import pandas as pd

def argparser():
parser = argparse.ArgumentParser(description="")
parser.add_argument("-maf", "--mafin", help="MAF input file", required=True)
parser.add_argument("-bed", "--bedout", help="BED input file", required=True)
parser.add_argument(
"--extra", help="Extra columns to keep (space separated list)", nargs="+", required=False, default=[]
)
return parser.parse_args()

def maf2bed(maf_file, bed_file, extra):
maf = pd.read_csv(maf_file, sep="\t", comment="#")
bed = maf[["Chromosome", "Start_Position", "End_Position"] + extra]
bed.to_csv(bed_file, sep="\t", index=False, header=False)

def main():
args = argparser()
maf2bed(maf_file=args.mafin, bed_file=args.bedout, extra=args.extra)

if **name** == "**main**":
main()

```

### 2. Make sure your script is in the right folder

Now, let's move it to the correct directory and make sure it is executable:

```bash
mv maf2bed.py /path/where/pipeline/is/bin/.
chmod +x /path/where/pipeline/is/bin/maf2bed.py
```

### 3. Create your custom module

Then, let's write our module. We will call the process
"CONVERT_MAF2BED" and add any tags or/and labels that
are appropriate (this is optional) and directives (via
conda and/or container) for
the definition of dependencies.

<blockquote style="border-left: 4px solid #F0AD4E; background-color: #FFF3CD; padding: 10px;">

<h4 style="margin-top: 0;">Some additional infos that might be of interest</h4>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on labels</span></summary>
A `label` will
annotate the processes with a reusable identifier of your
choice that can be used for configuring. E.g. we use the
`label` 'process_single', this looks as follows:

```

withLabel:process_single {
cpus = { check_max( 1 _ task.attempt, 'cpus' ) }
memory = { check_max( 1.GB _ task.attempt, 'memory') }
time = { check_max( 1.h \* task.attempt, 'time' ) }
}

```

</details>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on tags</span></summary>

A `tag` is simple a user provided identifier associated to
the task. In our process example, the input is a tuple
comprising a hash of metadata for the maf file called
`meta` and the path to the `maf` file. It may look
similar to: `[[id:'123', data_type:'maf'], /path/to/file/example.maf]`. Hence, when nextflow makes
the call and `$meta.id` is `123` name of the job
will be "CONVERT_MAF2BED(123)". If `meta` does not have
`id` in its hash, then this will be literally `null`.

</details>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on conda/container directives</span></summary>

The `conda` directive allows for the definition of the
process dependencies using the [Conda package manager](https://docs.conda.io/en/latest/). Nextflow automatically sets up an environment for the given package names listed by in the conda directive. For example:

```

process foo {
conda 'bwa=0.7.15'

'''
your_command --here
'''
}

```

Multiple packages can be specified separating them with a blank space e.g. `bwa=0.7.15 samtools=1.15.1`. The name of the channel from where a specific package needs to be downloaded can be specified using the usual Conda notation i.e. prefixing the package with the channel name as shown here `bioconda::bwa=0.7.15`

```

process foo {
conda 'bioconda::bwa=0.7.15 bioconda::samtools=1.15.1'

'''
your_bwa_cmd --here
your_samtools_cmd --here
'''
}

```

Similarly, we can apply the `container` directive to execute the process script in a [Docker](http://docker.io/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) container. When running Docker, it requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the local executor or the cluster nodes when the pipeline is deployed through a grid executor.

```

process foo {
conda 'bioconda::bwa=0.7.15 bioconda::samtools=1.15.1'
container 'dockerbox:tag'

'''
your_bwa_cmd --here
your_samtools_cmd --here
'''
}

```

Additionally, the `container` directive allows for a more sophisticated choice of container and if it Docker or Singularity depending on the users choice of container engine. This practice is quite common on official nf-core modules.

```

process foo {
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
'biocontainers/fastqc:0.11.9--0' }"

'''
your_fastqc_command --here
'''
}

```

</details>

</blockquote>

Since `maf2bed.py` is in the `bin` directory we can directory call it in the script block of our new module `CONVERT_MAF2BED`. You only have to be careful with how you call variables (some explanations on when to use `${variable}` vs. `$variable`):
A process may contain any of the following definition blocks: directives, inputs, outputs, when clause, and the process script. Here is how we write it:

```
process CONVERT_MAF2BED {
// HEADER
tag "$meta.id"
    label 'process_single'
    // DEPENDENCIES DIRECTIVES
    conda "anaconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
'quay.io/biocontainers/pandas:1.4.3' }"

// INPUT BLOCK
input:
tuple val(meta), path(maf)

// OUTPUT BLOCK
output:
tuple val(meta), path('\*.bed') , emit: bed
path "versions.yml" , emit: versions

// WHEN CLAUSE
when:
task.ext.when == null || task.ext.when

// SCRIPT BLOCK
script: // This script is bundled with the pipeline in bin
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"

"""
maf2bed.py --mafin $maf --bedout ${prefix}.bed
"""
}
```

More on nextflow's process components in the [docs](https://www.nextflow.io/docs/latest/process.html).

### Include your module in the workflow

In general, we will call out nextflow module `main.nf` and save it in the `modules` folder under another folder called `conver_maf2bed`. If you believe your custom script could be useful for others and it is potentially reusable or calling a tool that is not yet present in nf-core modules you can start the process of making it official adding a `meta.yml` [explained above](#adding-modules-to-a-pipeline). In the `meta.yml` The overall tree for the pipeline skeleton will look as follows:

```

pipeline/
├── bin/
│ └── maf2bed.py
├── modules/
│ ├── local/
│ │ └── convert_maf2bed/
│ │ ├── main.nf
│ │ └── meta.yml
│ └── nf-core/
├── config/
│ ├── base.config
│ └── modules.config
...

```

To use our custom module located in `./modules/local/convert_maf2bed` within our workflow, we use a module inclusions command as follows (this has to be done before we invoke our workflow):

```bash title="workflows/demotest.nf"
include { CONVERT_MAF2BED } from './modules/local/convert_maf2bed/main'
workflow {
input_data = [[id:123, data_type='maf'], /path/to/maf/example.maf]
CONVERT_MAF2BED(input_data)
}
```

:::tip{title="Exercise 6 - Adding a custom module"}
In the directory `exercise_6` you will find the custom script `print_hello.py`, which will be used for this and the next exercise.

1.  Create a local module that runs the `print_hello.py` script
2.  Add the module to your main workflow
3.  Run the pipeline
4.  Lint the pipeline
5.  Commit your changes
    <details>
    <summary>solution 1</summary>

    ```

    ```

      </details>

:::

<p class="text-center">
  <a href="/docs/contributing/nf_core_basic_training/add_nf_core_module/" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to Chapter 4
  </a>
  <a href="/docs/contributing/nf_core_basic_training/nf_schema/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 6 >
  </a>
</p>
