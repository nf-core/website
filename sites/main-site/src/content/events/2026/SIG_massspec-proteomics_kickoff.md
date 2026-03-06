---
title: "Mass Spec Proteomics SIG: Kickoff Meeting"
subtitle: "Initial meeting to discuss quantms integration, module development, and legacy pipelines."
type: "talk"
startDate: "2026-03-05"
startTime: "17:00+01:00"
endDate: "2026-03-05"
endTime: "18:00+01:00"
locations:
  - name: Virtual
---

## Overview

Following the recent introduction of the [Mass Spec Proteomics SIG](https://nf-co.re/blog/2025/sig-massspec-proteomics-introduction), we held our official kickoff meeting to discuss the current landscape of proteomics in nf-core, module development, and our roadmap for pipeline integration. It is incredibly encouraging to see the growing interest and enthusiasm for standardizing mass spectrometry and proteomics workflows within the nf-core community.

## Attendees

- [Phil Ewels](https://phil.ewels.co.uk/projects/nf-core/)
- [Martin Gordon](https://www.linkedin.com/in/martin-gordon-983a7a1b6/)
- [Dongze He](https://www.linkedin.com/in/dongzehe/)
- [enryh](https://orbit.dtu.dk/en/persons/henry-emanuel-webel/)
- SIG Lead

## Key Discussion Points

### 1. Collaboration and Integration with `quantms`

We spent a significant portion of the meeting discussing the excellent work done by the `bigbio` organization on the [`quantms`](https://github.com/bigbio/quantms) pipeline. While `quantms` is actively maintained in its own organization, we see a massive opportunity for mutual benefit by introducing their custom modules to the broader nf-core community.

**Current Progress:**

- We have successfully implemented the [`diann` module](https://github.com/nf-core/modules/tree/master/modules/nf-core/diann).
- The [`dia_proteomics_analysis` subworkflow](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/dia_proteomics_analysis) is currently available in nf-core.
- **Validation Success:** We confirmed that running the nf-core `thermorawfileparser` combined with the `dia_proteomics_analysis` subworkflow produces _exactly the same results_ as `quantms`. This includes identical DIANN configurations and step-by-step outputs.
- We are working on a comprehensive module and subworkflow collection for [fragpipe](https://fragpipe.nesvilab.org/) as well. This will be released in the upcoming weeks and will be a PR to nf-core/modules.

### 2. Expanding the Module Library

Our immediate goal is to transfer more modules from `quantms` into `nf-core/modules`. Key targets include:

- **OpenMS** modules ([https://github.com/bigbio/quantms/tree/master/modules/local/openms](https://github.com/bigbio/quantms/tree/master/modules/local/openms))
- **pmultiqc** module ([https://github.com/bigbio/quantms/tree/master/modules/local/pmultiqc](https://github.com/bigbio/quantms/tree/master/modules/local/pmultiqc))

**Objective:** By migrating these components to `nf-core/modules`, we aim to demonstrate to the `quantms` team that the shared nf-core ecosystem is highly robust. Relying on community-maintained modules can significantly reduce the maintenance burden for their pipeline developers, creating a win-win scenario.

### 3. The Future of Archived nf-core Proteomics Pipelines

We reviewed the status of several legacy nf-core proteomics pipelines that are currently archived:

- [`nf-core/diaproteomics`](https://github.com/nf-core/diaproteomics)
- [`nf-core/ddamsproteomics`](https://nf-co.re/ddamsproteomics/dev/docs/output/)
- [`nf-core/proteomicslfq`](https://github.com/nf-core/proteomicslfq)

**Next Steps:** The SIG will evaluate these pipelines to determine the community appetite and technical feasibility of reviving them versus directing efforts entirely toward modular components and new workflows.
