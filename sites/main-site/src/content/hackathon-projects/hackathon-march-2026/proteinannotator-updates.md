---
title: nf-core/proteinannotator updates
category: pipelines
slack: "https://nfcore.slack.com/archives/C084HP11VPF"
location: Wellcome Trust Campus
image: ""
image_alt: ""
leaders:
  vagkaratzas:
    name: Evangelos Karatzas
    slack: "https://nfcore.slack.com/team/U05LNHCFLCW"
  mberacochea:
    name: Martin Beracochea
    slack: "https://nfcore.slack.com/team/U04LW9A52MP"
---

This project aims to update the [nf-core/proteinannotator](https://nf-co.re/proteinannotator/dev) pipeline, by adding features described in the [issues](https://github.com/nf-core/proteinannotator/issues) page.

## Goal

Update the logic and add new features to the pipeline.

## Tasks

See below a couple of potential new features to add during the hackathon.
Depending on the participation and level of the attendees, more issues could be made available for the hackathon.

### Incorporate the NMPFamsDB in the domain annotation subworkflow

Issue [#77](https://github.com/nf-core/proteinannotator/issues/77)

[Novel Metagenome Protein Families (`NMPFamsDB`)](https://nmpfamsdb.org/) is a database which hosts novel metagenome protein clusters with no, or weak hits to Pfam or Reference genomes.
Use the HMM library from the resource, and search against it with `hmmer/hmmsearch` in the `domain_annotation` subworkflow, similarly to the existing `Pfam` and `FunFam` processes.

Level: _For anyone interested. Beginners welcome!_

### Add nf-test and meta.yml files for the local functional_annotation subworkflow

Issue [#78](https://github.com/nf-core/proteinannotator/issues/78)

The `domain_annotation` local subworkflow can be used as an example.

Level: _For anyone interested. Beginners welcome!_

### Incorporate the eggnogmapper software in the functional annotation subworkflow

Issue [#16](https://github.com/nf-core/proteinannotator/issues/16)

This should be done in two steps, similarly to interproscan;
first, download the required DBs, and second, run the eggnogmapper module.
The first task would be to update the existing `nf-core/modules` `EGGNOGMAPPER` module.
The next step would be to decide if it is better to create an nf-core module to download/preprocess the required databases, or re-use existing nf-core modules such as `ARIA2` and `UNTAR`.
There is already an existing PR ([#63](https://github.com/nf-core/proteinannotator/issues/63)) which can be salvaged for parts of the logic, or merge `dev` into that branch and continue from there.

Level: _For people with some experience in Nextflow/nf-core pipelines._
