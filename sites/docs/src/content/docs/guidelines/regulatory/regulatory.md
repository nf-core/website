---
title: Overview
subtitle: Overview on regulatory aspects for nf-core pipelines
weight: 1
---

# Introduction

nf-core aims to produce high-quality Nextflow pipelines that should make it easy to perform validation and adhere to regulatory requirements. The community as such is organized and discusses these requirements and sets these standards within the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory). Everyone is welcome to join in and bring in different points to it.

# Scope

## Developer & User guidelines

###  Versioning

### Code Contributions

* Protection of pipelines, everyone can contribute but releases are fixed / permanent
* Reviews
* Implementation of changes
  * Documentation (link to bug reports, suggestions, security updates)

### Documentation

* General documentation is available and complete
* Covers at least general aspects of the pipeline, e.g. how to run a basic example
* Covers the steps that are involved in the execution of a specifici subworkflow of the pipeline if there are multiple paths available within a pipeline

##  Testing

* nf-core tests, full-pipeline tests, unit tests with nf-test
* Process-level testing: test that every step in your workflow functions as expected
* Functional testing: run a full analysis test, based on a chosen set of test cases


## Maintenance

* Continous development
* Collection of bug reports, suggestions for new functionality
* User communication / communication guidelines

## Validation & integrative testing

* Integration testing: an analysis of the test in the production environment with real data
* Scope of this - what would this usually encompass?
* What is included in verification of performance specification of my test
system? You must compare the performance of the test system in your laboratory with those
  specifications established by the manufacturer. Define the following performance
  characteristics, on available test data:
  • Accuracy
  • Precision
  • Reportable range
  • Reference intervals/range
