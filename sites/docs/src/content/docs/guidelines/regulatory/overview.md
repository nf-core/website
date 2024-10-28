---
title: Overview
subtitle: Overview on regulatory aspects for nf-core pipelines
weight: 1
---

# Introduction

nf-core aims to produce high-quality Nextflow pipelines that should make it easy to perform validation and adhere to regulatory requirements. The community as such is organized and discusses these requirements and sets these standards within the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory). Everyone is welcome to join in and bring in different points to it.

# Scope

- Risk of the software application in the business process: What is the risk?
- (needs more work to be done)

## General Checklist

### Community

- History
- Size
- Governance model
- Users & Maintainers (number + experience of people)
- Issue tracking (number of issues and response time, closing time and ratio, classification / labeling) --> ask stats team
- Feature requests

### Versioning

- We enforce semantic versioning of pipeline releases in general
- Need a release - cannot work on a `dev` branch --> follow guidlines for release (link)
- Software tools inside a pipeline are fixed with respective docker/singularity containers that encapsulate (freeze) the tool used --> version is documented with pipeline release

### Code and software development process quality

- Protection of pipelines, everyone can contribute but releases are fixed / permanent
- How trustworthy is the pipeline?
- PR/code review and approval (change management)
- CI/CD (test automation, code coverage, versioning and releases)
- Patching and udpates (including frequency, monitoring of vulnerabilities and third party libraries), requirements management and technical documentation (traceability, reusability, granularity, updates)

- Reviews
- Implementation of changes
  - Documentation (link to bug reports, suggestions, security updates)

### Documentation

- General documentation is available and complete
- Covers at least general aspects of the pipeline, e.g. how to run a basic example
- Covers the steps that are involved in the execution of a specifici subworkflow of the pipeline if there are multiple paths available within a pipeline

## Testing

- nf-core tests, full-pipeline tests, unit tests with nf-test
- Process-level testing: test that every step in your workflow functions as expected
- Functional testing: run a full analysis test, based on a chosen set of test cases

## Maintenance

- Continous development
- Collection of bug reports, suggestions for new functionality
- User communication / communication guidelines
- Nf-core template updates create a new minor release at minimum --> not just a patch release
