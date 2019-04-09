---
title: Troubleshooting
subtitle: How to troubleshoot commong mistakes and issues
---


## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4. If you are running Single end data make sure to specify `--singleEnd`

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will.

## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Image cannot be build

Sometimes, `mksquashfs` cannot be found on the login node or workstation that you intend to use, thus the Singularity Image build fails unfortunately. See below code snippet that shows such a typical failure:
```bash
Caused by:
  Failed to pull singularity image
  command: singularity pull --name nfcore-rnaseq-1.3.img docker://nfcore/rnaseq:1.3 > /dev/null
  status : 255
  .....
    INFO:    Creating SIF file...
    FATAL:   Unable to pull docker://nfcore/rnaseq:1.3: While searching for mksquashfs: exec: "mksquashfs": executable file not found in $PATH
```
In such a case, please install `mksquashfs` or ask in your IT department that this package is installed for you. 

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us via the [Slack](https://nf-core-invite.herokuapp.com/) channel or by opening an issue in the respective pipeline repository on GitHub asking for help.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework [tools](https://github.com/nf-core/tools) then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
