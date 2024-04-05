---
title: Module Specifications
subtitle: Specifications for writing nf-core Nextflow DSL2 modules
---

## New module specifications

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### 1 General

#### 1.1 Required and optional input files

All mandatory and optional input files MUST be included in `input` channel definitions.

#### 1.2 Non-file mandatory command arguments

Non-file mandatory arguments or arguments needed to modify the command to make the module run with no error, SHOULD be provided as value channels (for example `lib_type` in [salmon/quant](https://github.com/nf-core/modules/blob/master/modules/nf-core/salmon/quant/main.nf)) - see 'Input/output options' below.

#### 1.3 Optional command arguments

All _non-mandatory_ command-line tool _non-file_ arguments MUST be provided as a string via the `$task.ext.args` variable.

- The value of `task.ext.args` is supplied from the `modules.config` file by assigning a closure that returns a string value to `ext.args`.
  The closure is necessary to update parameters supplied in a config with `-c`.

  ```groovy {2} title="<module>.nf"
  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  fastqc \\
      $args \\
        <...>
  """
  ```

  ```groovy {3-6} title="modules.config"
    process {
      withName: <module> {
          ext.args = { [                                                        // Assign a closure which returns a string
              '--quiet',
              params.fastqc_kmer_size ? "-k ${params.fastqc_kmer_size}" : ''    // Parameter dependent values can be provided like so
          ].join(' ') }                                                         // Join converts the list here to a string.
          ext.prefix = { "${meta.id}" }                                         // A closure can be used to access variables defined in the script
      }
    }
  ```

<details markdown="1">
  <summary>Rationale</summary>
  A disadvantage of passing arguments via ext.args is that it splits up how information is passed to a module, which can be difficult to understand where module inputs are defined.

The justification behind using the `ext.args` is to provide more flexibility to users.
As `ext.args` is derived from the configuration (e.g. `modules.config`), advanced users can overwrite the default `ext.args` and supply their own arguments to modify the behaviour of a module.
This can increase the capabilities of a pipeline beyond what the original developers intended.

Initially these were passed via the main workflow script using custom functions (e.g. `addParams`) and other additional nf-core custom methods, but this had a syntax overhead and other limitations that were found to be more difficult to use and understand by pipeline developers.
Therefore using the 'native' `ext` functionality provided by Nextflow was easier to understand, maintain and use.

Note that sample-specific parameters can still be provided to an instance of a process by storing these in `meta`, and providing these to the `ext.args` definition in `modules.config`.
A closure is used to make Nextflow evaluate the code in the code in the string.

```nextflow
ext.args = { "--id ${meta.id}" }
```

  </details>

#### 1.4 Use of multi-command piping

Software that can be piped together SHOULD be added to separate module files
unless there is a run-time, storage advantage in implementing in this way.

For example,
using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

```bash
bwa mem $args | samtools view $args2 -B -T ref.fasta
```

:::info
The addition of multi-tool modules to nf-core/modules adds increased burden on the nf-core
maintainers.
Where possible, if a multi-tool module is desired, it should be implemented as a local module in the nf-core pipeline.
If another nf-core pipeline also desires to use this module, a PR can be made adding it to nf-core/modules.
For guidelines regarding multi-tool modules, please search this page for the phrase `multi-tool`.
Existing local multi-tool modules can be searched for using the Github search box, searching across the nf-core org for terms such as `args2` `samtools` `collate` `fastq`.

```plaintext
org:nf-core args2 samtools collate fastq
```

Modules intended to batch process files by parallelizing repeated calls to a tool, for example with
`xargs` or `parallel`, also fall under the category of multi-tool modules.
Multi-tool modules
should chain tools in an explicit order given by the module name, e.g. `SAMTOOLS/COLLATEFASTQ`.
:::

#### 1.5 Each command must have an $args variable

Each tool in a multi-tool module MUST have an `$args` e.g.,

```bash
bwa mem $args | samtools view $args2 -B -T ref.fasta | samtools sort $args3
```

or

```bash
<tool> \\
   <subcommand> \\
   $args
gzip \\
    $args2
```

The numbering of each `$args` variable MUST correspond to the order of the tools, and MUST be documented in `meta.yml`.
E.g. in the first example, `bwa mem` is the first tool so is given `$args`, `samtools view` is the second tool so is `$args2`, etc.

#### 1.6 Types of meta fields

Modules MUST NOT use 'custom' hardcoded `meta` fields.
The only accepted 'standard' meta fields are `meta.id` or `meta.single_end`.
Proposals for other 'standard' fields for other disciplines must be discussed with the maintainers team.

<details markdown="1">
  <summary>Rationale</summary>
  Modules should be written to allow as much flexibility to pipeline developers as possible.

Hardcoding `meta` fields in a module will reduce the freedom of developers to use their own names for metadata, which would make more sense in that particular context.

As all non-mandatory arguments MUST go via `$args`, pipeline developers can insert such `meta` information into `$args` with whatever name they wish.

So, in the module code we DO NOT do:

```bash
my_command -r ${meta.strand} input.txt output.txt
```

... but rather, in `modules.conf`

```nextflow
ext.args = { "--r ${meta.<pipeline_authors_choice_of_name>}" }
```

... and then in the module code `main.nf`:

```bash
my_command $args input.txt output.txt
```

</details>

#### 1.7 Compression of input and output files

Where applicable, the usage and generation of compressed files SHOULD be enforced as input and output, respectively:

- `*.fastq.gz` and NOT `*.fastq`
- `*.bam` and NOT `*.sam`

If a tool does not support compressed input or output natively, we RECOMMEND passing the uncompressed data via unix pipes, such that it never gets written to disk, e.g.

```bash
gzip -cdf $input | tool | gzip > $output
```

The `-f` option makes `gzip` auto-detect if the input is compressed or not.

If a tool cannot read from STDIN, or has multiple input files, it is possible to use named pipes:

```bash
 mkfifo input1_uncompressed input2_uncompressed
 gzip -cdf $input1 > input1_uncompressed &
 gzip -cdf $input2 > input2_uncompressed &
 tool input1_uncompressed input2_uncompressed > $output
```

Only if a tool reads the input multiple times, it is required to uncompress the file before running the tool.

#### 1.8 Emission of versions

Where applicable, each module command MUST emit a file `versions.yml` containing the version number for each tool executed by the module, e.g.

```bash
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
END_VERSION
```

resulting in, for instance,

```yaml
'FASTQC':
  fastqc: 0.11.9
  samtools: 1.12
```

All reported versions MUST be without a leading `v` or similar (i.e. must start with a numeric character), or for unversioned software, a Git SHA commit id (40 character hexadecimal string).

<details class="mb-3">
<summary>Tips for extracting the version string</summary>

`sed{:bash}` is a powerful stream editor that can be used to manipulate the input text into the desired output.
Start by piping the output of the version command to `sed{:bash}` and try to select the line with the version number:

```bash
tool --version | sed '1!d'
```

- `sed '1!d'{:bash}` Extracts only line 1 of the output printed by `tools --version{:bash}`.
- The line to process can also be selected using a pattern instead of a number: `sed '/pattern/!d'{:bash}`, e.g. `sed '/version:/!d'{:bash}`.
- If the line extraction hasn't worked, then it's likely the version information is written to stderr, rather than stdout.
  In this case capture stderr using `|&{:bash}` which is shorthand for `2>&1 |{:bash}`.
- `sed 's/pattern/replacement/'{:bash}` can be used to remove parts of a string. `.` matches any character, `+` matches 1 or more times.
- You can separate `sed{:bash}` commands using `;`. Often the pattern : `sed 'filter line ; replace string'{:bash}` is enough to get the version number.
- It is not necessary to use `echo`, `head`, `tail`, or `grep`.
- Use `|| true` for tools that exit with a non-zero error code: `command --version || true{:bash}` or `command --version | sed ... || true{:bash}`.

</details>

We chose a [HEREDOC](https://tldp.org/LDP/abs/html/here-docs.html) over piping into the versions file line-by-line as we believe the latter makes it easy to accidentally overwrite the file.
Moreover, the exit status of the sub-shells evaluated in within the HEREDOC is ignored, ensuring that a tool's version command does no erroneously terminate the module.

If the software is unable to output a version number on the command-line then a variable called `VERSION` can be manually specified to provide this information e.g. [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf).

Please include the accompanying comments above the software packing directives and beside the version string.

```nextflow {4,15,21}
process TOOL {

...
// WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/tool:0.9.1--pl526hc9558a2_3' :
   'biocontainers/tool:0.9.1--pl526hc9558a2_3' }"

...

script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def VERSION = '0.9.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
"""
...

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool: $VERSION
END_VERSIONS
"""

}
```

If the HEREDOC cannot be used because the script is not bash, the `versions.yml` MUST be written directly e.g. [ascat module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf).

#### 1.9 Presence of when statement

The process definition MUST NOT change the `when` statement.
`when` conditions can instead be supplied using the `process.ext.when` directive in a configuration file.

```groovy
process {
    withName: 'FOO' {
        ext.when = !params.skip_module
    }
    withName: 'BAR' {
        ext.when = { meta.single_end }
    }
}
```

#### 1.10 Capturing STDOUT and STDERR

In some cases, STDOUT and STDERR may need to be saved to file, for example for reporting purposes.
Use the shell command `tee` to simultaneously capture and preserve the streams.
This allows for the streams to be captured by the job scheduler's stream logging capabilities and print them to screen when Nextflow encounters an error.

This also ensures that they are captured by Nextflow.

If information is only written to files, it could potentially be lost when the job scheduler gives up the job allocation.

```nextflow {7-8}
script:
"""
tool \\
  --input $input \\
  --threads $task.cpus \\
  --output_prefix $prefix \\
  2> >( tee ${prefix}.stderr.log >&2 ) \\
  | tee ${prefix}.stdout.log
"""
```

Similarly, if the tool captures STDOUT or STDERR to a file itself, it is best to send those to the corresponding streams as well.
Since a timeout may mean execution is aborted, it may make most sense to have background tasks do that.

```nextflow {3-4}
script:
"""
tail -F stored_stderr.log >&2 &
tail -F stored_stdout.log &
tool arguments
"""
```

#### 1.11 Capturing exit codes

Occasionally, some tools do not exit with the expected exit code 0 upon successful use of the tool.
In these cases one can use the `||` operator to run another useful command when the exit code is not 0 (for example, testing if a file is not size 0).

```nextflow {6}
script:
"""
tool \\
  --input $input \\
  --summary ${prefix}.summary \\
  || test -s ${prefix}.summary
"""
```

See the [Bash manual on file operators](https://tldp.org/LDP/abs/html/fto.html) for examples of properties of files which could be tested.

Alternate suggestions include using `grep -c` to search for a valid string match, or other tool which will appropriately error when the expected output is not successfully created.

### 2 Naming conventions

### 2.1 Name format of subworkflow files

The directory structure for the module name must be all lowercase, and without punctuation, e.g. [`modules/nf-core/bwa/mem/`](https://github.com/nf-core/modules/tree/master/modules/nf-core/bwa/mem/). The name of the software (i.e. `bwa`) and tool (i.e. `mem`) MUST be all one word.

Note that nf-core/tools will validate your suggested name.

#### 2.2 Name format of module processes

The process name in the module file MUST be all uppercase e.g. `process BWA_MEM {`. The name of the software (i.e. `BWA`) and tool (i.e. `MEM`) MUST be all one word separated by an underscore.

#### 2.3 Name format of module parameters

All parameter names MUST follow the `snake_case` convention.

#### 2.4 Name format of module functions

All function names MUST follow the `camelCase` convention.

#### 2.5 Name format of module channels

Channel names MUST follow `snake_case` convention and be all lower case.

#### 2.6 Command file output naming

Output file (and/or directory) names SHOULD just consist of only `${prefix}` and the file-format suffix (e.g. `${prefix}.fq.gz` or `${prefix}.bam`).

- This is primarily for re-usability so that other developers have complete flexibility to name their output files however they wish when using the same module.
- As a result of using this syntax, if the module has the same named inputs and outputs then you can add a line in the `script` section like below (another example [here](https://github.com/nf-core/modules/blob/e20e57f90b6787ac9a010a980cf6ea98bd990046/modules/lima/main.nf#L37)) which will raise an error asking the developer to change the `args.prefix` variable to rename the output files so they don't clash.

  ```nextflow
  script:
  if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
  ```

### 3 Input/output options

#### 3.1 Required `path` channel inputs

Input channel `path` declarations MUST be defined for all _possible_ input files (i.e. both required and optional files).

- Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel (e.g. [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
- Other generic auxiliary files used across different input files (e.g. common reference sequences) MAY be defined using a dedicated input channel (e.g. [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

#### 3.2 Required `val` channel inputs

Input channel `val` declarations SHOULD be defined for all mandatory non-file inputs that are essential for the functioning of the tool (e.g. parameters, flags etc).

- Mandatory non-file inputs are options that the tool MUST have to be able to be run.
- These non-file inputs are typically booleans or strings, and must be documented as such in the corresponding entry in the `meta.yaml`.
- Options, flags, parameters that are _not_ required by the tool to function should NOT be included - rather these can be passed via `ext.args`.

 <details markdown="1">
 <summary>Rationale</summary>
 It was decided by a [vote](https://nfcore.slack.com/archives/C043UU89KKQ/p1677581560661679) amongst interested parties within the 2023 Maintainers group on 2023-02-28 to allow non-file mandatory input channels.

The reasoning behind this was that it is important to have documented (using the existing display on the website) the bare minimum information required for a module to run.
It also allows module code to consume parameter values without parsing them out of the `ext.args` string and reduces possible risks of entire breakage of modules with future [expected config changes](https://github.com/nextflow-io/nextflow/issues/2723) at a Nextflow level.

Downsides to this approach are readability (now multiple places must be checked on how to modify a module execution - modules.conf `ext.args`, the module invocation in pipeline code etc.), and reduced user freedom.
However it was felt that it was more important for stability in and 'installation' and 'execution' of modules was preferred (e.g. for tools that require position arguments etc.)

</details>

<details markdown="1">
<summary>Inputs particular cases</summary>
 When one and only one of multiple argument are required:

- If they all are string argument : use 1 argument that will be equal to the string

  e.g. Parameter model of [glimpse2 chunk](https://nf-co.re/modules/glimpse2_chunk)

- If some are files put them all in one channel and test if only one is present

  e.g. Grouping output parameters of [glimpse2 concordance](https://nf-co.re/modules/glimpse2_concordance)

  `if (((file1 ? 1:0) + (val1 ? 1:0) + (val2 ? 1:0)) != 1) error "One and only one argument required"`
  </details>

#### 3.3 Output channel emissions

Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

#### 3.4 Optional inputs

Optional inputs are not currently supported by Nextflow.
However, passing an empty list (`[]`) instead of a file as a module parameter can be used to work around this issue.

For example, having a module (`MY_MODULE`) that can take a `cram` channel and an optional `fasta` channel as input, can be used in the following ways:

```nextflow
MY_MODULE(cram, [])     // fasta is optional, the module will run without the fasta present
MY_MODULE(cram, fasta)  // execution of the module will need an element in the fasta channel
```

#### 3.5 Optional outputs

Optional outputs SHOULD be marked as optional:

```nextflow
tuple val(meta), path('*.tab'), emit: tab,  optional: true
```

#### 3.6 One output channel per output file type

Each output file SHOULD be emitted in its own channel (and no more than one), along with the `meta` map if provided ( the exception is the versions.yml ).

### 4 Documentation

#### 4.1 Module documentation is required

Each module MUST have a `meta.yaml` in the same directory as the `main.nf` of the module itself.

#### 4.2 Number of keywords

Keywords SHOULD be sufficient to make the module findable through research domain, data types, and tool function keywords

- Keywords MUST NOT just be solely of the (sub)tool name

:::info
For multi-tool modules, please add the keyword `multi-tool`, as well as all the (sub)tools involved.
:::

#### 4.3 Keyword formatting

Keywords MUST be all lower case

#### 4.4 Documenting of all tools

The tools section MUST list every tool used in the module. For example

```yml
tools:
  - bowtie2: <....>
  - samtools: <....>
```

#### 4.5. Documentation of args of each piped or multiple command

The tools section MUST have a `args_id:` field for every tool in the module that describes which `$args` (`$args2`, `$args3`) variable is used for that specific module. A single tool module will only have `args_id: "$args"`.

```yml
tools:
  - bowtie2:
      <...>
      args_id: "$args"
  - samtools:
      <...>
      args_id: "$args2"
```

#### 4.6 Required channel documentation

Input and Output sections of the `meta.yaml` SHOULD only have entries of input and output channels

#### 4.7 Documentation of tuples

Input and output tuples MUST be split into separate entries

- i.e., `meta` should be a separate entry to the `file` it is associated with

#### 4.8 Input and output channel types

Input/output types MUST only be of the following categories: `map`, `file`, `directory`, `string`, `boolean`, `integer`, `float`, `boolean`, `list`

#### 4.9 Correspondence of input/outputs entries to channels

Input/output entries MUST match a corresponding channel in the module itself

- There should be a one-to-one relationship between the module and the `meta.yaml`
- Input/output entries MUST NOT combine multiple output channels

#### 4.10 Useful input/output descriptions

Input/output descriptions SHOULD be descriptive of the contents of file

- i.e., not just 'A TSV file'

#### 4.11 Input/output glob pattern

Input/output patterns (if present) MUST follow a [Java glob pattern](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob)

#### 4.12 Indication of input channel requirement

Input entries should be marked as Mandatory or Optional

### 5 Module parameters

#### 5.1 Module input and outputs

A module file SHOULD only define input and output files as command-line parameters to be executed within the process.

#### 5.2 Use of parameters within modules

All `params` within the module MUST only be initialised and used in the local context of the module.
In other words, named `params` defined in the parent workflow MUST NOT be assumed to be passed to the module to allow developers to call their parameters whatever they want.
In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

#### 5.3 Specification of multiple-threads or cores

If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

#### 5.4. Evaluation of parameter within a module

Any parameters that need to be evaluated in the context of a particular sample e.g. single-end/paired-end data MUST also be defined within the process.

### 6 Resource requirements

#### 6.1 Use of labels in modules

An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config#L29-L46) e.g. `process_single`, `process_low`, `process_medium` or `process_high`.

#### 6.2. Source of multiple threads or cores value

If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

If the tool does not support multi-threading, consider `process_single` unless large amounts of RAM are required.

#### 6.3 Specifying multiple threads for piped commands

If a module contains _multiple_ tools that supports multi-threading (e.g. [piping output into a samtools command](https://github.com/nf-core/modules/blob/c4cc1db284faba9fc4896f64bddf7703cedc7430/modules/nf-core/bowtie2/align/main.nf#L47-L54)), you can assign CPUs per tool.

- Note that [`task.cpus`] is supplied unchanged when a process uses multiple cores
- If one tool is multi-threaded and another uses a single thread, you can specify directly in the command itself e.g. with [`${task.cpus}`](https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/sampe/main.nf#L34)

### 7 Software requirements

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/).
Where possible we will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

#### 7.1 Use of container directives

Software requirements SHOULD be declared within the module file using the Nextflow `container` directive.
For single-tool BioContainers, the `nf-core modules create` command will automatically fetch and fill-in the appropriate Conda / Docker / Singularity definitions by parsing the information provided in the first part of the module name:

```nextflow
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    'biocontainers/fastqc:0.11.9--0' }"
```

#### 7.2 Use of conda directive

If the software is available on Conda it MUST also be defined in an `environment.yml` file alongside the `main.nf` of the module, and is passed to the Nextflow `conda` directive within `main.nf`.

Using `bioconda::bwa=0.7.17` as an example, software MUST be pinned to the channel (i.e. `bioconda`) and version (i.e. `0.7.17`).

Conda packages MUST not be pinned to a build because they can vary on different platforms.

#### 7.3 Re-use of multi-tool containers

If required, multi-tool containers may also be available on BioContainers e.g. [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40).
You can install and use the [`galaxy-tool-util`](https://anaconda.org/bioconda/galaxy-tool-util) package to search for both single- and multi-tool containers available in Conda, Docker and Singularity format.
E.g. to search for Docker (hosted on Quay.io) and Singularity multi-tool containers with both `bowtie` and `samtools` installed you can use the following command:

```console
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep "mulled"
```

:::note
Build information for all tools within a multi-tool container can be obtained in the `/usr/local/conda-meta/history` file within the container.
:::

#### 7.4 Creation of new multi-tool containers

It is also possible for a new multi-tool container to be built and added to BioContainers by submitting a pull request on their [`multi-package-containers`](https://github.com/BioContainers/multi-package-containers) repository.

- Fork the [multi-package-containers repository](https://github.com/BioContainers/multi-package-containers)
- Make a change to the `hash.tsv` file in the `combinations` directory see [here](https://github.com/aunderwo/multi-package-containers/blob/master/combinations/hash.tsv#L124) for an example where `pysam=0.16.0.1,biopython=1.78` was added.
- Commit the code and then make a pull request to the original repo, for [example](https://github.com/BioContainers/multi-package-containers/pull/1661)
- Once the PR has been accepted a container will get built and you can find it using a search tool in the `galaxy-tool-util conda` package

  ```console
  mulled-search --destination quay singularity conda  --search pysam biopython  | grep "mulled"
  quay         mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  docker pull quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
  singularity  mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  wget https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
  ```

- You can copy and paste the `mulled-*` path into the relevant Docker and Singularity lines in the Nextflow `process` definition of your module
- To confirm that this is correct, spin up a temporary Docker container

  ```console
  docker run --rm -it quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0  /bin/sh
  ```

  And in the command prompt type

  ```console
  $ grep specs /usr/local/conda-meta/history
  # update specs: ['biopython=1.78', 'pysam=0.16.0.1']
  ```

  The packages should reflect those added to the multi-package-containers repo `hash.tsv` file

- If the multi-tool container already exists and you want to obtain the `mulled-*` path, you can use (this)[https://midnighter.github.io/mulled] helper tool.

#### 7.5 Software not on Bioconda

If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. We will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).

### 8 Misc

#### 8.1 General module code formatting

All code MUST be aligned to follow the '[Harshil Alignment™️](#what-is-the-harshil-alignment)' format.
