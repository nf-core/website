---
title: Stuck on revision
subtitle: How to troubleshoot common mistakes and issues
weight: 50
---

If you get a warning like the following:

```bash
Project nf-core/<pipeline> is currently stuck on revision: dev -- you need to explicitly specify a revision with the option -r in order to use it
```

This is a Nextflow error, with less-commonly seen Git 'terminology'. What this means is that you have multiple versions of the pipeline pulled (e.g. 2.0.0, 2.1.0, 2.1.1, dev etc.), and it is not sure which one to use. Therefore, with every `nextflow run nf-core/<PIPELINE>` command you should always indicate which version with `-r`.

For example:

```bash
nextflow run nf-core/<pipeline> -r 2.1.0 --input '/<path>/<to>/data/*_{R1,R2}_*.fq.gz' <...>
```

Specifying the version of the run you are using is highly recommended, as it helps in full reproducibility. In the sense that if you explicitly record the whole command _with_ the version for your publication or internal reports, then anyone who wants to check your work can use the exact version you used (including all internal tools).

You can see more information on the Nextflow documentation [here](https://www.nextflow.io/docs/latest/sharing.html#handling-revisions).

:::note
This warning previous had various typos of `stickied` or `sticked`, so depending on which version of Nextflow you're using you may see that mentioned.
:::
