---
title: nf-core/proteinfold Extension
category: pipelines
location: Skopje (Loka Office)
leaders:
  - JelPej
---

Extending the existing [nf-core/proteinfold](https://nf-co.re/proteinfold/1.0.0/)
pipeline to support protein complex co-folding and structural validation metrics.

## Goal

Add three new tools to nf-core/proteinfold to enable co-folding of protein complexes
and evaluation of prediction quality — making the pipeline more powerful for
drug discovery and structural biology research.

## Tasks

1. **Add Boltz-2 co-folding module**

   Write a new nf-core module for `boltz/predict`.
   Boltz-2 supports protein-protein, antibody-antigen and protein-ligand co-folding,
   and predicts both structure and binding affinity.


   Level: Intermediate — some Nextflow/nf-core experience recommended.

2. **Add US-align structural comparison module**

   Write a new nf-core module
   for `usalign/align`. US-align measures RMSD and TM-score between a predicted
   structure and a known reference — quantifying how accurate the prediction is.

 
   - Source: https://github.com/pylelab/USalign
 

   Level: Intermediate — some Nextflow/nf-core experience recommended.

3. **Add DockQ v2 interface scoring module**

   Write  new nf-core module for `dockq/score`.
   DockQ specifically scores the quality of the binding interface between two proteins 

 
   - Source: https://github.com/bjornwallner/DockQ

   Level: Intermediate — some Nextflow/nf-core experience recommended.

4. **Wire new modules into the pipeline**

   Connect Boltz-2 → US-align → DockQ into the existing proteinfold workflow
   and validate with the built-in test profile.

   Level: Intermediate — Nextflow pipeline development experience recommended.

