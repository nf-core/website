---
title: Overview
subtitle: Overview on regulatory aspects for nf-core pipelines
weight: 1
---

# Introduction

nf-core aims to produce high-quality Nextflow pipelines that should make it easy to perform validation and adhere to regulatory requirements. The community as such is organized and discusses these requirements and sets these standards within the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory). Everyone is welcome to join in and bring in different points to it.

# Scope

* Risk of the software application in the business process: What is the risk?
* License
*

## Developer & User guidelines


### Community

* History
* Size
* Governance model
* Users & Maintainers (number + experience of people)
* Issue tracking (number resopnse time, closing time and ratio, classification / labeling) --> ask stats team
* Feature requests


###  Versioning

###  Code and software development process quality

* Protection of pipelines, everyone can contribute but releases are fixed / permanent
* How trustworthy is the pipeline?
* PR/code review and approval (change management)
* CI/CD (test automation, code coverage, versioning and releases)
* Patching and udpates (including frequency, monitoring of vulnerabilities and third party libraries), requirements management and technical documentation (traceability, reusability, granulatirty, updates)

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
* How is potential risk derived from aspects listed above and what mitigation measures are implemented?
