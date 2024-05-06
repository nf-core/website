---
title: Harshil Alignment™️
subtitle: How to format your Nextflow code to make nf-core core members happy
---

The Harshil Alignment™️ format is the whitespace-happy code style that was introduced by a certain core member to get on everyone's nerves, but then make subsequently develop Stockholm Syndrome so that no-one in nf-core else now can look at Nextflow code without it.

The Harshil Alignment™️ format involves ensuring that common punctuation across multiple lines in a group are placed in the same location as each other.

There are many places where the format can be applied - it's not just code, it can also applies to comment formatting - however common examples are as follows:

### Curly Bracket Example

❌ Bad

```groovy
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
```

✅ Good

```groovy
include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
```

### Equals Example

❌ Bad

```groovy
stats = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
```

✅ Good

```groovy
stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats    ] ]
flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
```

### Comma Example

❌ Bad

```groovy
tuple val(meta), path("*.bam"), emit: bam, optional:true
tuple val(meta), path("*.log"), emit: log
tuple val(meta), path("*fastq.gz"), emit: fastq, optional:true
path  "versions.yml", emit: versions
```

✅ Good

```groovy
tuple val(meta), path("*.bam")    , emit: bam     , optional:true
tuple val(meta), path("*.log")    , emit: log
tuple val(meta), path("*fastq.gz"), emit: fastq   , optional:true
path  "versions.yml"              , emit: versions
```

### Colon Example (Comments)

```groovy
take:
print_version        // boolean: print version
dump_parameters      // boolean: dump parameters
outdir               // path: base directory used to publish pipeline results
check_conda_channels // boolean: check conda channels
```

✅ Good

```groovy
take:
print_version        // boolean: print version
dump_parameters      // boolean: dump parameters
outdir               //    path: base directory used to publish pipeline results
check_conda_channels // boolean: check conda channels
```
