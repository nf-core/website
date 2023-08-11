---
title: Module Writing Guidelines
subtitle: Guidelines on writing a sanger-tol compliant module
---

# NF-Core Module Guidelines

Generally, modules should conform to the [NF-Core Standards](https://nf-co.re/docs/contributing/modules#new-module-guidelines-and-pr-review-checklist). However, there are a small number of cases where this is not possible.

# Sanger-tol additions
## Example 1: The Super-module

In [TreeVal](https://www.github.com/sanger-tol/treeval), there is the [`cram_filter_align_bwamem2_fixmate_sort.nf`](https://github.com/sanger-tol/treeval/blob/dev/modules/local/cram_filter_align_bwamem2_fixmate_sort.nf) local module, which contains 3 tools across 5 commands.

```bash
input:
tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix)

...

script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
cram_filter -n ${from}-${to} ${cramfile} - | \\
    samtools fastq -F0xB00 -nt - | \\
    bwa-mem2 mem -p ${bwaprefix} -t${task.cpus} -5SPCp -H'${rglines}' - | \\
    samtools fixmate -mpu - - | \\
    samtools sort --write-index -l1 -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -
"""
```

The original implementation of this module took a cram file and whilst reading, split the stream into 10,000 container segments (a cram container is analogous to 1 read), these stream segments are then processed upon (by the other 4 commands) before being merged into 1 mapped and sorted bam file for the primary assembly.

Nextflow cannot manipulate a data stream passing between two modules. This required us to create a module to pre compute the 10,000 container regions of interest in the cram file (in the form of a csv) and pass these as arguments to the cram_et al_ module. Whilst not as performant as the original implementation, this is much more performant (in terms of compute resources and IO impact) than splitting the cram file into n (( total number of container / 10,000 ) * no. of cram files) number of files before further manipulation with the next 4 commands. This means that the TreeVal implementation is the best case scenario, as shown below.

<img src="../../public_html/assets/img/developer-images/cram-et-al.png" alt="A comparison of the three different cram et al module implementations, first the original (and fastest), second the TreeVal and finally a wholly NF-Core implementation" width="600" height="500">

### Reasons for:

- Where there is a significant reduction in compute resource and IO when compared to a subworkflow with the same function.

### Reasons against:

- Where there is no performance increase.
- Where an intermediary file, inside the command pipe, is also an output file.
- More complex testing (requires more resources as, in this case, there are 5 tools running concurrently).
- Mulled containers are more cryptic as to the packages contained.
