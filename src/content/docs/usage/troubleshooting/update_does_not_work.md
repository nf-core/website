---
title: Pipeline update doesn't work
subtitle: How to troubleshoot common mistakes and issues
---

### My pipeline update doesn't seem to do anything

To download new version of a pipeline, you can use the following.

```bash
nextflow pull nf-core/<pipeline> -latest
```

To download a previous version of a pipeline, you can instead use the following, replacing `<version>` to the corresponding version (v3.4 for example)

```bash
nextflow pull nf-core/<pipeline> -r `<version>`
```

However, in very rare cases, minor fixes to a version will be pushed out without a version number bump. This can confuse Nextflow slightly, as it thinks you already have the 'broken' version from your original pipeline download.

> This _shouldn't_ happen with stable versions and normally only happens on `dev` branches.

If when running the pipeline you don't see any changes in the fixed version when running it, you can try removing your Nextflow's nf-core pipeline cache typically stored in your home directory with:

```bash
rm -r ~/.nextflow/assets/nf-core/<pipeline>
rm -r ~/.config/nfcore/nf-core/modules
```

And re-pull the pipeline with the command above. This will install a fresh version of the version with the fixes.
