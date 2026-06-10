---
title: nf-core/cellpainting
category: pipelines
slack: https://nfcore.slack.com/archives/C08P26ZAZM4
image: "/assets/images/events/2026/hackathon-boston/cellpainting_meme.png"
image_alt: "The 'X all the Y' meme stick figure shouting wildly with the caption 'PAINT ALL THE CELLS'."
leaders:
  kenibrewer:
    name: Ken Brewer
    slack: https://nfcore.slack.com/team/U03MWA6LMV1
---

This project focuses on [nf-core/cellpainting](https://github.com/nf-core/cellpainting), a Nextflow pipeline for scalable, reproducible image-based profiling of [Cell Painting](https://doi.org/10.1038/nprot.2016.105) assays. The pipeline takes high-content microscopy images, applies illumination correction, extracts morphological features with [CellProfiler](https://cellprofiler.org/), and converts the results to analysis-ready [Parquet](https://parquet.apache.org/) files with [CytoTable](https://github.com/cytomining/CytoTable).

The pipeline already has a working core (illumination correction → assay-development QC → CellProfiler analysis → CytoTable conversion → MultiQC), but it is pre-1.0 and there is a long list of high-value features and polish work to do before a first release. The Boston hackathon is a great chance to broaden the contributor base and knock out a batch of well-scoped issues across the stack.

---

## Goal

Make meaningful progress toward a 1.0 release of nf-core/cellpainting by closing a batch of [open issues](https://github.com/nf-core/cellpainting/issues) — adding new analysis modules from the cytomining ecosystem, tightening the test suite, and improving downstream usability of the pipeline outputs.

---

## What participants will do

Each contributor will:

1. Pick an issue from the [nf-core/cellpainting issue tracker](https://github.com/nf-core/cellpainting/issues) and assign themselves.
2. Discuss the approach with the project lead and other contributors at the table.
3. Implement the change on a feature branch (new module, test improvement, or documentation update).
4. Run nf-test and lint locally with the `test` profile.
5. Open a Pull Request for review.

---

## Suggested tasks

### Good first issue

- [#40 — Switch CytoTable grouping from per-site to per-plate](https://github.com/nf-core/cellpainting/issues/40). Self-contained Nextflow refactor: regroup the CytoTable invocation by `(batch, plate)` and update the test snapshot. Touches one module, one workflow file, and the docs.

### Test infrastructure

- [#39 — Add a `test_full` profile that runs a full plate from the minimal test dataset](https://github.com/nf-core/cellpainting/issues/39). Generate a full-plate samplesheet for `BR00117035`, host it on the `cellpainting` branch of `nf-core/test-datasets`, and wire up the AWS megatest profile.
- [#38 — Replace coarse md5 snapshots with intelligent assertions for CellProfiler outputs](https://github.com/nf-core/cellpainting/issues/38). Move from blanket file-ignore lists to nf-test path/line assertions that allow non-deterministic fields (run timestamps, work-dir paths, PNG metadata) without losing coverage.
- [#15 — Add SQLite Cell Painting Gallery example data](https://github.com/nf-core/cellpainting/issues/15). Subset the `BR00117035.sqlite` from the Cell Painting Gallery and contribute it to `nf-core/test-datasets`.

### New modules — cytomining downstream stack

- [#9 — `pycytominer_annotate`](https://github.com/nf-core/cellpainting/issues/9), [#10 — `pycytominer_aggregate`](https://github.com/nf-core/cellpainting/issues/10), [#11 — `pycytominer_normalize`](https://github.com/nf-core/cellpainting/issues/11), [#12 — `pycytominer_featureselect`](https://github.com/nf-core/cellpainting/issues/12). Add the [pycytominer](https://github.com/cytomining/pycytominer) operations as nf-core modules and wire them into the analysis subworkflow downstream of CytoTable.
- [#8 — `cosmicqc`](https://github.com/nf-core/cellpainting/issues/8). Add [coSMicQC](https://github.com/WayScience/coSMicQC) as a single-cell QC step.

### New modules — imaging

- [#2 — `cellpose`](https://github.com/nf-core/cellpainting/issues/2). Integrate [Cellpose](https://www.cellpose.org/) segmentation as an alternative to CellProfiler primary/secondary objects, feeding label matrices into a re-configured analysis pipeline.
- [#5 — `fiji_stitchsegmentedimages`](https://github.com/nf-core/cellpainting/issues/5). Stitch per-site assay-development overlays into a pseudo-plate view for at-a-glance visual QC.

---

## Who should join?

This project is a good fit for people interested in:

- Nextflow / nf-core pipeline development (any experience level)
- Image-based profiling, high-content imaging, or Cell Painting assays
- The cytomining ecosystem (CellProfiler, CytoTable, pycytominer, coSMicQC)
- Test engineering with nf-test
- Bringing a new pipeline to its first release

You do **not** need a biology background — most of the open issues are pure pipeline-engineering tasks.

---

## Recommended preparation

- Basic familiarity with Git and GitHub (forks, branches, pull requests)
- Some exposure to Nextflow and nf-core ([Hello nextflow](https://training.nextflow.io/latest/hello_nextflow/) and [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/) are great primers)
- Working Docker or Singularity install for running `-profile test,docker`
- Optional reading: the [Cell Painting protocol](https://doi.org/10.1038/nprot.2016.105) and the [JUMP Cell Painting Consortium](https://jump-cellpainting.broadinstitute.org/) overview for context on the data
