---
title: User and developer aspects - regulatory 
subtitle: Developer and user specific documentation for regulatory aspects in nf-core
weight: 1
---


## Scope 

This document intends to provide further context for maintainers and developers. 
## Validation & integrative testing

* Scope :
   This part of the document provides information on which features should be included with each pipeline, which would facilitate regulatory approvals:
* Features to include:
    Pipeline level:
  * Integration testing: an analysis of the test in the production environment with real data
    - Compare the performance of the test system in your dataset with those
    specifications defined by the user. This includes the following performance
    characteristics:
        • Accuracy
        • Precision
    if applicable    • Reportable range
    if applicable    • Reference intervals/range (normal values) for the laboratory’s patient population
* [if applicable] Controls to be included to unit tests:
        • Positive control
        • Negative control
        • Additional controls (for example PCR reagent controls, amplification control gene, calibration curve,... )
* Set of expected results for all controls.
* Set assay acceptance criteria
* Set rejection criteria.
   - Add automated quality checks that will be stopping points for the pipeline, if fail.

* Store results to an automated report / stats file
* Automate risk management based on results stored in the stats file
  
