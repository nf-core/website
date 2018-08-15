# How does Nextflow work?
The general outline of a Nextflow script is as follows, and it will be described in detail in the following sections:

```bash
#!/usr/bin nextflow

params.first = < first.input.parameter >
params.sec = < second.input.parameter >

process < name > {

  [ directives ]

  input:
   < process inputs >

  output:
   < process outputs >

  script:
  < user script to be executed >

}

```
## Processes and Channels
Each script usually has two or more [**processes**](https://www.nextflow.io/docs/latest/process.html), which are steps of the pipeline that work independent from each other.
For example could the first process produce library indexes and the second process could map reads to the indexed library.


Information can be moved from one processes of a pipeline to another via [**channels**](https://www.nextflow.io/docs/latest/channel.html).
They are pieces of information defined in the *input* and *output* section of each process and can be files, parameters, lists, etc.


Interestingly, the order of processes within a nextflow script is not important as the workflow is defined by input and output.
A process will be initiated once all the elements of the channels defined in the input are available.


## Scripts
Each nextflow process ***must end*** with a **script** block, which defines the code to be executed by the process.
It can call on the variables from the input section and will release the output to the output section.
A script is initiated as follows:

```
script:
"""
<user script>
"""
```

The `script:` header can be omitted if using triple quote blocks, however we recommend always using it for clarity.

Note that a script block delimited by a `"` character support variable substitutions, while blocks delimited by `'` do not.

It is possible to use any type of programming language for execution, simply by adding a leading shebang statement to the top of the `script` block (eg. `!#/usr/bin/env perl`).
If no shebang is specified, the command is performed using bash.


# Executors
In the Nextflow framework architecture, the [**executor**](https://www.nextflow.io/docs/latest/executor.html) is the component that determines the system where a pipeline process is run and supervises its execution.
The executor provides an abstraction between the pipeline processes and the underlying execution system.
This allows you to write the pipeline functional logic independently from the actual processing platform.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager or the cloud by simply changing the executor definition in the Nextflow configuration file.


# Nextflow tutorial
The following steps will guide you through the setup of a typical Nextflow RNAseq pipeline with practical examples to follow.

The [`nextflow-io/crg-course-nov16` GitHub repository](https://github.com/nextflow-io/crg-course-nov16) contains the tutorial material for the *Parallel distributed computational workflows
with Nextflow and Docker containers* course.

Clone this repository with the following command in your terminal:

```bash
git clone https://github.com/nextflow-io/crg-course-nov16.git
cd crg-course-nov16
```
Make sure Nextflow and Docker is installed on your computer.
You can find instructions about this in the [usage documentation](/usage_docs).

## Nextflow hands-on

During this tutorial you will implement a proof of concept of a RNA-Seq pipeline which:

1. Indexes a genome file.
2. Maps read pairs against the genome.
3. Performs quantification.

### Step 1 - Command line parameters

The script `rna-ex1.nf` defines the pipeline input parameters that can be defined on the command line but does not contain any processes yet.
Such parameters follow the convention `params.<name>` in the `run rna-ex1.nf` script file (see line 5-7 for examples).

This is how the parameter definition looks in the script:
```groovy
params.genome = "$baseDir/data/ggal/genome.fa"
```

Because [Nextflow is based on Groovy](https://www.nextflow.io/docs/latest/script.html#language-basics), you will find groovy functions like `println()` in the scripts (see lines `12` to `16`).

Run `rna-ex1.nf` with the default parameters by using the following command:

```bash
nextflow run rna-ex1.nf
```

To specify a different [input parameter](https://www.nextflow.io/docs/latest/config.html#scope-params), you can use the parameter in the command prefixed with a double dash (`--`). For example:

```bash
nextflow run rna-ex1.nf --genome "this/and/that"
```
In this example the path leading to the genome file `"$baseDir/data/ggal/genome.fa"` in the script will be overwritten with `"this/and/that"` when you run the command.

> Note: Parameter arguments must be enclosed with quotes if they contain spaces or a file glob pattern

### Step 2 - Build genome index

`rna-ex2.nf`contains the first [process](https://www.nextflow.io/docs/latest/process.html) called `buildIndex`. It takes the genome file as
input and creates the genome index by using the `bowtie2-build` tool.

You may recall, that `params.genome` contains only the path to a file.
In order to access the file we need to evoke the [`file` method](https://www.nextflow.io/docs/latest/script.html?highlight=file#files-and-i-o) (see line `22` of `rna-ex2.nf`):

```groovy
genome_file = file(params.genome)
```

This is then used in the `buildIndex` process:
```groovy
input:
file genome from genome_file
```

In this example the object `genome_file` is used as a [channel](https://www.nextflow.io/docs/latest/channel.html) for the input and connects the `params.genome` to the variable `genome` inside the process.
As such, the variable `genome` can be used in the process to call the information from the channel `genome_file`.
This variable can only be used within this process, _except_ if it is declared in the output.

The following command will fail because [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is not installed in the test environment.
Try to run it by using the command:

```bash
nextflow run rna-ex2.nf
```

Add the command line option `-with-docker` to launch the execution through a Docker container as shown below:

```bash
nextflow run rna-ex2.nf -with-docker
```

This time it works because it uses the Docker container `nextflow/rnatoy:1.3` defined in [the
`nextflow.config` file](https://www.nextflow.io/docs/latest/config.html#configuration-file).
The `config.file`  in Nextflow is a simple text file containing a set of properties defined using the syntax: `name = value`

In order to avoid having to specify the option `-with-docker` every time, you can add the following line in the `nextflow.config` file:

```groovy
docker.enabled = true
```

You will find that the process also contains an `output` block as well as a `script` block.
In the `script` block, bowtie2 is executed in a bash environment using the variable  `genome` from the input block as input file.
The `bowtie2` output will be passed to the `output` block via the variable `genome_index`.
The generated bowtie2 output file is named as defined after `file` in the output block, in this case `genome.index<extension>` (see lines 33 to 34).

When a process is initiated it is reported on the terminal with a specific hash identifier like `[22/7548fa]` which [identifies the unique process execution](https://www.nextflow.io/docs/latest/tracing.html?highlight=task#tracing-visualisation). This code is also the prefix of the directories where each process is executed.
You can inspect the files produced in the work subdirectory generated automatically in the Nextflow project directory.

> Note: Nextflow has built-in default variables.
> In this script we see the example of [`${task.cpus}`](https://www.nextflow.io/docs/latest/tracing.html?highlight=task#trace-report), which defines the number of cpus requested for a task execution.
 > As these variables are already defined, they do not need to be initiated in the script, but they can be overwritten in either the config file or in the process.

### Step 3 - Collect read files by pairs

This step shows how to match paired-end FastQ files into pairs, so they can be mapped by [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml).
The script `rna-ex3.nf` adds the creation of the `read_pairs` channel.
The [operator `.fromFilePairs`](https://www.nextflow.io/docs/latest/operator.html) emits tuples containing three elements: the pair ID, the first read-pair file and the second read-pair file.

Edit the script `rna-ex3.nf` and add the following statement as the last line:

```groovy
read_pairs.println()
```

Save it and execute it with the following command:

```bash
nextflow run rna-ex3.nf
```

Try it again specifying different read files by using a glob pattern:

```bash
nextflow run rna-ex3.nf --reads 'data/ggal/reads/*_{1,2}.fq'
```

It shows how read files matching the pattern specified are grouped in pairs having
the same prefix.

> NB: The quotes used with `--reads` here are needed to stop your bash command line from automatically expanding the glob path before it gets to nextflow. Without them, the read pairing will not work.



### Step 4 - Map sequence reads

The script `rna-ex4.nf` adds the `mapping` process.
Note how it declares three inputs:
the genome fasta file, the genome index file produced by the `buildIndex` process and
the read pairs.
Also note that the last [input is defined utilizing `set`](https://www.nextflow.io/docs/latest/process.html?highlight=set#inputs) allowing it to define three different elements: a pair ID, the first read file and the second read file.
Reminiscent to `set` in the input block, in the output block, the [`set` qualifier](https://www.nextflow.io/docs/latest/process.html?highlight=set#output-set-of-values) allows to send multiple values into a single channel.
This allows to group together the results of multiple executions of the same process.

Execute `rna-ex4.nf` by using the following command:

```bash
nextflow run rna-ex4.nf -resume
```

The `-resume` option skips the execution of any step that has been processed in a previous
execution.

Try to execute it with more read files as shown below:

```bash
nextflow run rna-ex4.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```


### Step 5 - Perform reads quantification

In script `rna-ex5.nf` the `makeTranscript` process is added in order to quantify the reads.
It takes the annotation file defined in the parameters at the beginning of the nextflow script and the *bam* files produced by *TopHat* in the mapping process.
Then it quantifies the reads using [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/).
Cufflinks in it’s default settings produces an output file names `transcripts.gtf`.
In order to keep the sample information, the script renames the files adding the pair ID by using the bash command `mv`.

> Note: The script block contains two lines of bash code. This is possible because we use triple `"""`.

You can run `rna-ex5.nf` by using the following command:

```bash
nextflow run rna-ex5.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```

### Step 6 - Define the pipeline output

In script `rna-ex6.nf` the [`publishDir` directive](https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir) is added in order to produce the pipeline output to a folder of your choice.
`mode : 'copy'` allows us to copy the files into the chosen directory while keeping the original files in the corresponding `work` subfolder.

Run the example by using the following command:

```bash
nextflow run rna-ex6.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```

Then you will find the quantification files in the folder `results`.

Modify the `rna-ex6.nf` script by adding the following line at the beginning of the file:

```groovy
params.outdir = 'results'
```

Then, look for the `publishDir` directive in the `makeTranscript` process, and
replace the `'results'` string with the `params.outdir` parameter.
By doing so we can change the destination folder for the final output from the command line.

Run `rna-ex6.nf ` again with the following command:

```bash
nextflow run rna-ex6.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq' --outdir my_transcripts
```

You will find the transcripts produced by the pipeline in the `my_transcripts` folder.


### Step 7 - Handle completion event

This step shows how to execute an action when the pipeline completes the execution.
In script `rna-ex7.nf` the [`tag` directive](https://www.nextflow.io/docs/latest/process.html?highlight=tag#tag) allows you to associate each process executions with a custom label, making easier to identify them in the log file or in the trace execution report.

> Note: Nextflow processes define the execution of *asynchronous* tasks. They are not
executed in the order that they are written in the pipeline script as it would happen in a
common *iterative* programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message
when the script completes.

Try to run it by using the following command:

```bash
nextflow run rna-ex7.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```

**If you want to expand your Nextflow skills in writing pipelines and managing Docker containers, we recommend to follow [this course](https://github.com/nextflow-io/rmghc-2018).**
