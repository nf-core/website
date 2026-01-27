---
title: General
subtitle: General guidelines for nf-core module development
markdownPlugin: addNumbersToHeadings
shortTitle: General
weight: 1
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Required and optional input files

All mandatory and optional input files MUST be included in `input` channel definitions.

## Non-file mandatory command arguments

Non-file mandatory arguments required for the module to run without error SHOULD be provided as value channels (e.g., `lib_type` in [salmon/quant](https://github.com/nf-core/modules/blob/master/modules/nf-core/salmon/quant/main.nf)).

## Optional command arguments

All _non-mandatory_ command-line tool _non-file_ arguments MUST be supported in the module via the `$task.ext.args` variable.

The `$args` variable MUST be placed in the module's tool command to allow optional and/or dynamic variables to be specified by a user or a developer.
The contents of this variable is specified by a pipeline developer or user in with the `ext.args` process variable in a `modules.config` or other Nextflow config file using a closure.
Config file are present in pipelines and subworkflows, not for each individual module, except for configs for tests.
The `$ext.args` variable can also be used in module tests to increase the test coverage of a tool's functionality through a companion `nextflow.config` alongside the test files themselves.

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
            params.fastqc_kmer_size ? "-k ${params.fastqc_kmer_size}" : ''    // Parameter-dependent values can be provided in this way
        ].join(' ') }                                                         // Join converts the list to a string.
        ext.prefix = { "${meta.id}" }                                         // A closure can be used to access variables defined in the script
    }
  }
```

:::info{title="Rationale" collapse}
Passing arguments via ext.args splits how information is passed to a module, making it harder to understand where module inputs are defined.

Using `ext.args` provides more flexibility to users.
As `ext.args` is derived from the configuration (e.g.,, `modules.config`), advanced users can overwrite the default `ext.args` and supply their own arguments to modify the behaviour of a module.
This can increase the capabilities of a pipeline beyond what the original developers intended.

Initially, these were passed via the main workflow script using custom functions (e.g., `addParams`) and other additional nf-core custom methods, but this syntax overhead and other limitations were more difficult for pipeline developers to use and understand.
Using the 'native' `ext` functionality provided by Nextflow was easier to understand, maintain, and use.

Sample-specific parameters can still be provided to an instance of a process by storing these in `meta`, and providing these to the `ext.args` definition in `modules.config`.
A closure is used to make Nextflow evaluate the code in the string:

```groovy
ext.args = { "--id ${meta.id}" }
```

:::

## Use of multi-command piping

Software that can be piped together SHOULD be added to separate module files unless this provides run-time or storage advantages.

For example, using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

```bash
bwa mem $args | samtools view $args2 -B -T ref.fasta
```

:::info
Multi-tool modules in `nf-core/modules` increase the burden on nf-core maintainers.
Where possible, implement multi-tool modules as local modules in the nf-core pipeline.
If another nf-core pipeline needs to use this module, make a PR to add it to nf-core/modules.
For guidelines regarding multi-tool modules, search this page for the phrase `multi-tool`.
Search for existing local multi-tool modules using the GitHub search box across the nf-core org for terms such as `args2` `samtools` `collate` `fastq`.

```plaintext
org:nf-core args2 samtools collate fastq
```

Modules intended to batch process files by parallelizing repeated calls to a tool, for example with
`xargs` or `parallel`, also fall under the category of multi-tool modules.
Multi-tool modules should chain tools in an explicit order given by the module name.
For example, `SAMTOOLS/COLLATEFASTQ`.
:::

## Each command must have an $args variable

Each tool in a multi-tool module MUST have an `$args`:

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
For example, in the first example, `bwa mem` is the first tool so is given `$args`, `samtools view` is the second tool so is `$args2`, etc.

## Types of meta fields

'Custom' hardcoded `meta` fields MUST NOT be used in modules.
Do not refer to them within the module as expected input, nor generate new fields as output.
The only accepted 'standard' meta fields are `meta.id` or `meta.single_end`.
Discuss proposals for other 'standard' fields for other disciplines with the maintainers team on slack under the [#modules channel](https://nfcore.slack.com/archives/CJRH30T6V).

:::info{title="Rationale" collapse}
Write modules to allow as much flexibility to pipeline developers as possible.

Hardcoding `meta` fields as input and output reduces the freedom of developers to use their own metadata names in their specific context.

As all non-mandatory arguments MUST go via `$args`, pipeline developers can insert such `meta` information into `$args` with whatever name they wish.

In the module code DO NOT:

```bash title="main.nf"
my_command -r ${meta.strand} input.txt output.txt
```

... but rather:

```groovy title="modules.conf"
ext.args = { "--r ${meta.<pipeline_authors_choice_of_name>}" }
```

and then in the module code:

```bash title="main.nf"
my_command $args input.txt output.txt
```

:::

## Compression of input and output files

Where applicable, compressed files SHOULD be used as input and output:

- `*.fastq.gz` and NOT `*.fastq`
- `*.bam` and NOT `*.sam`

If a tool does not support compressed input or output natively, nf-core RECOMMENDS passing the uncompressed data via unix pipes so that it never gets written to disk, for example:

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

Only if a tool reads the input multiple times, uncompress the file before running the tool.

## Emission of versions

The topic output qualifier in Nextflow collects outputs from multiple processes across a pipeline.
Use this feature in nf-core modules to collect version information from all tools without complex channel mixing logic.
See the [fastqc module](https://github.com/nf-core/modules/blob/0c47e4193ddde2c5edbc206b5420cbcbee5c9797/modules/nf-core/fastqc/main.nf#L16) as an example.

:::warning
For modules that use the template process directive, they will currently continue to depend on the old approach with `versions.yml`.
The only difference is that they should also use the topic output qualifier to send the `versions.yml` file to the versions topic.
:::

:::tip{title="Tips for extracting the version string" collapse}

`sed{:bash}` is a powerful stream editor that can be used to manipulate the input text into the desired output.
Start by piping the output of the version command to `sed{:bash}` and try to select the line with the version number:

```bash
tool --version | sed '1!d'
```

- `sed '1!d'{:bash}` Extracts only line 1 of the output printed by `tools --version{:bash}`.
- The line to process can also be selected using a pattern instead of a number: `sed '/pattern/!d'{:bash}`. For example, `sed '/version:/!d'{:bash}`.
- If the line extraction hasn't worked, then it's likely the version information is written to stderr, rather than stdout.
  In this case capture stderr using `|&{:bash}` which is shorthand for `2>&1 |{:bash}`.
- `sed 's/pattern/replacement/'{:bash}` can be used to remove parts of a string. `.` matches any character, `+` matches 1 or more times.
- You can separate `sed{:bash}` commands using `;`. Often the pattern: `sed 'filter line ; replace string'{:bash}` is enough to get the version number.
- It is not necessary to use `echo`, `head`, `tail`, or `grep`.
- Use `|| true` for tools that exit with a non-zero error code: `command --version || true{:bash}` or `command --version | sed ... || true{:bash}`.

:::

:::note
For not yet converted modules, you will see a different approach for collecting versions.
Even though the approach is deprecated, it is shown below for reference.
:::

Where applicable, each module command MUST emit a file `versions.yml` containing the version number for each tool executed by the module, for example:

```bash
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
END_VERSIONS
```

resulting in, for instance,

```yaml
"FASTQC":
  fastqc: 0.11.9
  samtools: 1.12
```

All reported versions MUST be without a leading `v` or similar (that is, must start with a numeric character), or for unversioned software, a Git SHA commit id (40 character hexadecimal string).

A [HEREDOC](https://tldp.org/LDP/abs/html/here-docs.html) is used over piping into the versions file line-by-line to avoid accidentally overwriting the file.
The exit status of sub-shells evaluated within the HEREDOC is ignored, ensuring that a tool's version command does not erroneously terminate the module.

If the software is unable to output a version number on the command-line, manually specify a variable called `VERSION` to provide this information.
For example, [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf).

Include the accompanying comments above the software packing directives and beside the version string.

```groovy {4,15,21}
process TOOL {

...
// WARN: Version information not provided by tool on CLI. Update version string below when bumping container versions.
conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/tool:0.9.1--pl526hc9558a2_3' :
   'biocontainers/tool:0.9.1--pl526hc9558a2_3' }"

...

script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def VERSION = '0.9.1' // WARN: Version information not provided by tool on CLI. Update this string when bumping container versions.
"""
...

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool: $VERSION
END_VERSIONS
"""

}
```

If the HEREDOC cannot be used because the script is not bash, write the `versions.yml` directly.
For example, [ascat module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf).

## Presence of when statement

The `when` statement MUST NOT be changed in the process definition.
Supply `when` conditions using the `process.ext.when` directive in a configuration file instead.

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

## Capturing STDOUT and STDERR

STDOUT and STDERR may need to be saved to a file for reporting or debugging.
The `tee` shell command captures and preserves these streams simultaneously.

This setup ensures that the job scheduler can capture stream logs while also printing them to the screen if Nextflow encounters an error.
This is especially useful when using `process.scratch` (which executes the process in a temporary folder), as logs might otherwise be lost on error.
The stream output is preserved in the process's `.command.log` and `.command.err` files.
If information is only written to files, it could be lost if the job scheduler reclaims the job allocation.

```groovy {7-8}
script:
"""
tool \\
  --input $input \\
  --threads $task.cpus \\
  --output_prefix $prefix \\
  2>| >( tee ${prefix}.stderr.log >&2 ) \\
  | tee ${prefix}.stdout.log
"""
```

:::tip{title="Reason for forced stream redirect" collapse}

nf-core sets the `-C` (`noclobber`) flag for each shell process, which prevents redirection from overwriting existing files.
Since some shells also treat this stream redirection as an error, we use the forced redirection `2>|` instead of `2>`.

:::

Similarly, if the tool itself captures STDOUT or STDERR to a file, it's best to redirect those to the corresponding streams as well.
For instance, if a timeout aborts execution, it's often more reliable to have background tasks handle this redirection.

```groovy {3-4}
script:
"""
tail -F stored_stderr.log >&2 &
tail -F stored_stdout.log &
tool arguments
wait
"""
```

## Capturing exit codes

Some tools do not exit with the expected exit code of 0 upon successful execution.
In these cases, use the `||` operator to run another useful command when the exit code is not 0 (for example, testing if a file is not size 0).

```groovy {6}
script:
"""
tool \\
  --input $input \\
  --summary ${prefix}.summary \\
  || test -s ${prefix}.summary
"""
```

See the [Bash manual on file operators](https://tldp.org/LDP/abs/html/fto.html) for examples of file properties you can test.

Alternative suggestions:

- Use `grep -c` to search for a valid string match
- Use another tool that errors when the expected output is not created

## Script inclusion

Module templates separate scientific logic from workflow-specific logic, improving code clarity and maintainability.
If a module's `script:` block contains a script rather than command invocations, regardless of the language (for example, Bash, R, Python), and the content is more than a readable length (as a rule of thumb, approximately 20 lines), provide it through a [Nextflow module template](https://www.nextflow.io/docs/latest/module.html#module-templates).

:::note
We recommend use of Nextflow templates as they are the most portable method of separating custom script content and execution across all execution contexts.
:::

:::note
Where script content in a module becomes particularly extensive, we strongly encourage packaging and hosting the code externally and provisioning via Conda/Docker as a standalone tool(kit).
:::

### Inline script code

If the script content remains at a readable length, the code MAY be embedded directly in the module without a dedicated template file.
However, embedded scripts should still follow the guidance for templates.

### Module template location

Place the template in a directory called `templates/` in the same directory as the module `main.nf`.

Name the template file after the module itself with a language-appropriate file suffix.
For example, the `deseq2/differential` nf-core module will use the `deseq2_differential.R`.

Refer to the template file within the module using the template function:

    ```nextflow
    script:
    template 'deseq2_differential.R'
    ```

See [`deseq2/differential`](https://github.com/nf-core/modules/blob/master/modules/nf-core/deseq2/differential/main.nf#L47) for an example of a template in an nf-core pipeline.

The resulting structure would look like this.

    ```tree
    deseq2
    └── differential
        ├── environment.yml
        ├── main.nf
        ├── meta.yml
        ├── templates
        │   └── deqseq2_differential.R
        └── tests
            ├── main.nf.test
            ├── main.nf.test.snap
            └── tags.yml
    ```

### Template or inline script-code contents

:::warning
Be aware that in any script template, Nextflow variables need to be escaped in the same way as in a standard bash `script:` block.
:::

The script template file or inline script code (used when at a readable length) MUST generate a `versions.yml` file using language-appropriate methods that contains versions of the base language and all relevant libraries and packages.

The generated `versions.yml` MUST have the same structure as a standard nf-core module `versions.yml`.

See the [`deseq2/differential` module](https://github.com/nf-core/modules/blob/4c2d06a5e79abf08ba7f04c58e39c7dad75f094d/modules/nf-core/deseq2/differential/templates/deseq_de.R#L509-L534) for an example using R.

### Stubs in templated modules

A templated module MUST have a stub block in the same way as any other module.
For example, use `touch` to generate empty files and versions.
See [`deseq2/differential` module](https://github.com/nf-core/modules/blob/4c2d06a5e79abf08ba7f04c58e39c7dad75f094d/modules/nf-core/deseq2/differential/main.nf#L34-L49) for an example in an nf-core module.

An inline command MAY be used to call the version for libraries for the `versions.yml` in this case.
For an R example see [deseq2/differential](https://github.com/nf-core/modules/blob/4c2d06a5e79abf08ba7f04c58e39c7dad75f094d/modules/nf-core/deseq2/differential/main.nf#L47).

## Stubs

### Stub block must exist

[A stub block](https://www.nextflow.io/docs/latest/process.html#stub) MUST exist for all modules.
This is a block of code that replaces the `script` command when the option `-stub` is set.
This enables quick testing of the workflow logic, as a "dry-run".

### Stub block prefix and versions

Include the same variables (for example, `prefix`) and HEREDOC code in the stub block as in the main script block.

### Stub files for all output channels

Include the creation of at least one file for every output channel (both mandatory and optional) in the stub block, generated with touch, for example:

```groovy
stub:
"""
touch ${prefix}.txt
"""
```

Ideally, the stub block should reproduce as much as possible the number of files and filename structure of the files expected as output.

### Stub gzip files must use echo and pipe

Use the syntax in the following example for stub files output as gzip compressed:

```bash
echo "" | gzip > ${prefix}.txt.gz
```

:::info{title="Rationale" collapse}
Touching a file with the file name ending in `.gz` will break nf-test's Gzip file parser, as the file is not actually gzipped and cannot be read.

Generate a valid gzipped file for nf-test to accept it during tests.
:::
