---
title: User and developer aspects - regulatory 
subtitle: Developer and user specific documentation for regulatory aspects in nf-core
weight: 1
---


## Scope 

This document intends to provide further context for maintainers and developers. 
## Validation & integrative testing

* Integration testing: an analysis of the test in the production environment with real data
* Scope :
    Pipeline level:
    - Compare the performance of the test system in your dataset with those
    specifications defined by the user. This includes the following performance
    characteristics:
        • Accuracy
        • Precision
    if applicable    • Reportable range
    if applicable    • Reference intervals/range (normal values) for the laboratory’s patient population
* Add automated quality checks that will be stopping points for the pipeline, if fail.
* Store results to an automated report / stats file
* Automate rish management based on results stored in the stats file
  
