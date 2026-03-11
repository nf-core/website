---
title: QuantMS pipeline (DIANN, nf-core updates, documentation, benchmarking)
category: pipelines
location: "ZS Copenhagen"
slack: "https://nfcore.slack.com/archives/C02Q3FL29PD"
leaders:
  enryh:
    name: Henry Webel
    slack: "https://nfcore.slack.com/archives/U07FM4513GD"
    github: "https://github.com/enryh"
---

# QuantMS pipeline

The goal is to make a state of affairs. Compatibility with nf-core, preparing
local for nf-core/modules (which can still be patched to be locally updated). We
as the organizing team are interested in running DIA-NN through the workflow, which
can be done either via QuantMS or the nf-core/dia_proteomics_analysis` subworkflow.

Additionally exploring the new tools `nf-docs` and `nf-metro` for documentation
and deployment instructions can be done during the hackathon.

Last, if people want to rather explore datasets, benchmarking of microbial datasets
can be done.

Good to prepare

- local up-to-date docker, nextflow and nf-core tools installation
- nf-core tutorial for modules and subworkflows

## DIA-NN updates

- document how to use a DIA-NN image using a private registry
- DIA-NN updates: [#663](https://github.com/bigbio/quantms/issues/663)
- try to use `ext.args` see [PR660](https://github.com/bigbio/quantms/pull/660)

### Subworkflows

Process names are not aligned, but I mapped them one to one in order.

- DIANN subworkflow in `nf-core/modules` repo at
  `subworkflows/nf-core/dia_proteomics_analysis`

  ```groovy
  include { QUANTMSUTILS_DIANNCFG          } from '../../../modules/nf-core/quantmsutils/dianncfg/main'
  include { QUANTMSUTILS_MZMLSTATISTICS    } from '../../../modules/nf-core/quantmsutils/mzmlstatistics/main'
  include { QUANTMSUTILS_DIANN2MZTAB       } from '../../../modules/nf-core/quantmsutils/diann2mztab/main'

  include { DIANN as DIANN_INSILICOLIBRARYGENERATION } from '../../../modules/nf-core/diann/main'
  include { DIANN as DIANN_PRELIMINARYANALYSIS } from '../../../modules/nf-core/diann/main'
  include { DIANN as DIANN_ASSEMBLEEMPIRICALLIBRARY } from '../../../modules/nf-core/diann/main'
  include { DIANN as DIANN_INDIVIDUALANALYSIS } from '../../../modules/nf-core/diann/main'
  include { DIANN as DIANN_FINALQUANTIFICATION } from '../../../modules/nf-core/diann/main'
  ```

- DIANN under local modules in `bigbio/quantms`. Process names are not aligned, but I
  mapped them one to one in order. So files could be compared to the ones in
  `nf-core/modules` repo.

  ```groovy
  include { GENERATE_CFG                } from '../modules/local/diann/generate_cfg/main'
  include { MSSTATS_LFQ                 } from '../modules/local/msstats/msstats_lfq/main'
  include { CONVERT_RESULTS             } from '../modules/local/diann/convert_results/main'

  include { INSILICO_LIBRARY_GENERATION } from '../modules/local/diann/insilico_library_generation/main'
  include { PRELIMINARY_ANALYSIS        } from '../modules/local/diann/preliminary_analysis/main'
  include { ASSEMBLE_EMPIRICAL_LIBRARY  } from '../modules/local/diann/assemble_empirical_library/main'
  include { INDIVIDUAL_ANALYSIS         } from '../modules/local/diann/individual_analysis/main'
  include { FINAL_QUANTIFICATION        } from '../modules/local/diann/final_quantification/main'
  ```

Compare and take inspiration by Jonathan Mannings way to write modules and subworkflows?

## Updates to nf-core

To get familiar with nf-core templates and requirements, one could try to move
some tools for the use of others to `nf-core/modules` repo. Any in

- `subworkflows/local`
- `modules/local`
- `modules/bigbio`

One could use and update modules which have a local version, but are maintained by others
in `nf-core/modules` repo. For example:

- Update ThermoRawFileParser (C#) to use `nf-core/modules/thermorawfileparser` version
  instead of `modules/bigbio/thermorawfileparser`

## Updates to nextflow language

In the [bigbio/nf-modules](https://github.com/bigbio/nf-modules) repo

- read and parse `args` and `prefix,` see
  [#23](https://github.com/bigbio/nf-modules/issues/23)
- fix issues shown by the language server, see
  [#24](https://github.com/bigbio/nf-modules/issues/24)

Compare with `nf-core/modules` repo, see
[here](https://github.com/bigbio/nf-modules/tree/main/modules/bigbio)

There are currently only 3 modules in
[bigbio/nf-modules](https://github.com/bigbio/nf-modules) repo, so further could be
added.

W.r.t to [#23](https://github.com/bigbio/nf-modules/issues/23) and
[#24](https://github.com/bigbio/nf-modules/issues/24)
we need to check:

- `thermorawfileparser`
- `onsite`

### Exercise: Add a module to nf-core/modules

- if the process is based on a python or conda package, wave allows easy
  containerization
- nf-tests need to be added

Useful hints.

- [create a template](https://nf-co.re/docs/nf-core-tools/modules/create).
- [check module specs](https://nf-co.re/docs/guidelines/components/modules)
- [test-data](https://nf-co.re/docs/guidelines/components/test_data)

List of candidates (tbc)

- pmultiqc (Python)
- msstats (R)

## nf-core lint

`.nf-core.yaml` file deactivates some things for linting. check what and how.

Run

```bash
# in quantms repo
nf-core pipelines lint -d .
```

## nf-docs

Add or use [ewels/nf-docs](https://github.com/ewels/nf-docs)

## nf-metro

Add a new metro-map based on a configuration file:
[pinin4fjords/nf-metro](https://pinin4fjords.github.io/nf-metro/latest/guide/)

- maybe add to deployment instructions (manual updates or actions)

## Comparisons using PRIDE: DIA datasets

Run experiments, compare outputs to results provided on PRIDE. Familiarize with
running quantms. Should be supplemented with inhouse data, which is now all DIA
on Bruker experiments.

- [PXD054415](https://www.ebi.ac.uk/pride/archive/projects/PXD054415) - comparing
  DDA and DIA on metaproteomics dataset with known compositions
  - could use a subset of samples
  - [x] SDRF
- [PXD049262](https://www.ebi.ac.uk/pride/archive/projects/PXD049262)
  - growth experiment
  - photosynthetic metabolism of purple sulphur bacteria Halorhodospira halophila
  - cultivated with various sulphur compounds
  - [ ] SDRF

### Included benchmark dataset in quantms

Mentioned as an example for
[DIA](https://docs.quantms.org/en/latest/benchmarks.html#data-independent-acquisition-dia)

- [article](https://www.sciencedirect.com/science/article/pii/S2352340922000415)
- [data](https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=69c2b1bd22cd4933887b4b4846da1bd7)

## Performance benchmarking

- running quantms on a single machine, single VM, on Azure batch, on HPCs with apptainers:
  - runtime, costs, etc.

## DIANN docker files

- update docker container to latest version of DIA-NN, test with apptainer:
  [bigbio/quantms-containers](https://github.com/bigbio/quantms-containers)
