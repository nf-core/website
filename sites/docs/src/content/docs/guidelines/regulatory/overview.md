---
title: Overview
subtitle: Overview on regulatory aspects for nf-core pipelines
weight: 1
---

# Introduction

nf-core aims to produce high-quality Nextflow pipelines that should make it easy to perform validation and adhere to regulatory requirements.
The community as such is organized and discusses these requirements within the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory).
Everyone is welcome to join in and bring in different points to it.

## Document Version

1.0.0 draft

## Scope

Pipeline validation is widely seen as something critical for use cases where outputs and interpretation of results matter more strictly, e.g. where regulatory authorities impose certain quality requirements to be met.
To assess this, the usual approach is to perform a risk based validation.

Risk based validation is considering everything around the development, implementation and integration of analysis pipelines as a potential risk or threat in terms of a misuse or malfunction of a pipeline.
This then has to be mitigated using appropriate measures.
An example would be that a pipeline per se has the risk of failing execution, which is a risk of it not producing desired outcome and which can be mitigated using appropriate functional tests using [nf-test](https://www.nf-test.com/) for example.

> [!WARNING]
> While nf-core can provide users with guidelines, information and help to validate pipelines, we will not be able to provide you with a full validation report that you can simply take "off the shelf" and use for your regulatory needs.

    The report that nf-core will be able to create for you [soon](https://github.com/nf-core/tools/issues/3258) will however contain a lot of the basic information required for running a full validation.

> We are working on a proof of concept validation for one nf-core pipeline (rnaseq) to showcase what needs to be done and where are potential gaps within the nf-core guidelines, processes or tooling that we can then hopefully address.

"Unless specifically exempted in a classification regulation, any medical device software product developed after June 1, 1997, regardless of its device class, is subject to applicable design control provisions." [FDA and IEEE definitions](https://www.complianceonline.com/resources/software-verification-and-validation-requirements-for-medical-device-and-implementation-strategies.html)

### What types of validation "exist"

Regulatory authorities play a pivotal role in the pharmaceutical industry, particularly in the context of drug approval.
These entities, such as the U.S. Food and Drug Administration (FDA) consensus standards, the European Medicines Agency (EMA)/ European Commission, and others, are responsible for ensuring the safety, efficacy, and quality of drugs and medical tests before they reach the market.
Each of them has its own set of standards that need to be followed, depending on the intended use of the test. We differentiate between standards that apply on a software level versus on a more infrastructure & computer systems validation level.
The latter (CSV) is out of scope for this guidance and should be dealt with depending on your infrastructure or IT provider.

Potentially applicable guidelines for bioinformatics pipelines:

- [SaMD general](https://www.greenlight.guru/blog/samd-software-as-a-medical-device)
- [FDA SaMD](https://www.fda.gov/medical-devices/digital-health-center-excellence/software-medical-device-samd)
- [CE mark registration](https://europa.eu/youreurope/business/product-requirements/labels-markings/ce-marking/index_en.htm)
- [FDA LDT](https://www.fda.gov/medical-devices/in-vitro-diagnostics/laboratory-developed-tests)
- [CLIA validation](https://www.cms.gov/regulations-and-guidance/legislation/clia/downloads/6064bk.pdf)
- [Medical device registration](https://health.ec.europa.eu/system/files/2021-10/mdcg_2021-24_en_0.pdf)

Computerized systems validation (out of scope):

- [CSV validation](https://en.wikipedia.org/wiki/Computerized_system_validation)

Because of the number of current standards, there are initiatives to harmonize them, see [Food and Drug Administration (FDA) consensus standards](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/appropriate-use-voluntary-consensus-standards-premarket-submissions-medical-devices) and [European Commission (EC) Harmonized Standards](https://single-market-economy.ec.europa.eu/single-market/european-standards/harmonised-standards_en) for references on this.

## Guidelines for validating nf-core pipelines

These points are individual points you should consider when validating an nf-core pipeline. A simpler version of this as a checklist is provided on a separate document page without the full explanations.

### Community metrics

From a risk-based perspective, open-source communities present both benefits and risks in validating analysis pipelines.
Consider the following metrics and criteria to judge the quality and risks associated with a pipeline:

- History of the pipeline
  - When was this started?
  - How many users started this?
  - How many people have contributed to this until now?
- Size
  - How complex is the pipeline? Does it involve many steps or very few steps? Are the building blocks (modules, subworkflows local or have they been reviewed / being worked upon within the larger community, e.g. stem from nf-core/modules)
- Governance model
  - We should reference here the nf-core governance webpage, which should also be versioned (see ticket on versioning of the website)
- Licensing
  - All nf-core pipelines utilize the MIT licence, which allows for public as well as commercial usage of said pipelines
  - Pipelines that do not adhere to this, have extra documentation that highlights which parts require specific software licences
- Users & Maintainers
  - How many active users have been using / cloning the pipeline in the last year?
  - How many people have contributed overall / in the last year to the pipeline?
  - How trustworthy are the main contributors? Have they been active a lot (set a threshold for this) in the community?
- Issue tracking
  - How many open issues do we have right now?
  - How many issues were opened overall, how many were closed?
  - Can we assess the quality of issues, are they ranked around severity?
  - How long does it take to close issues, how many are closed by the maintainers?
  - How is feature requesting handled within the pipeline maintainers?

### General requirements

- Define the functional requirements of your pipeline upfront
- Maintain a comprehensive list of any outside dependencies (tools, references, utilized public information, ...)
- Verify infrastructure requirements for consistent and timely pipeline execution. Note: Infrastructure validation itself falls under Computerized Systems Validation (CSV) and is out of scope for this document.

### Versioning

nf-core pipelines enforce semantic versioning for [pipeline releases](https://nf-co.re/docs/guidelines/pipelines/requirements/semantic_versioning).

Semantic versioning provides [a simple set of rules and requirements that dictate how version numbers are assigned and incremented](https://semver.org/). Version numbers have a meaning about how the underlying code has been modified from one version to another. To ensure quality and compliance, it is advised to conduct automated testing at appropriate levels aligned with the type of release. As a summary:

- Patches (x.y.Z | x > 0) introduce backward compatible bug fixes.
  Testing for these changes should focus on ensuring that existing functionality remains unaffected.
- Minor versions (x.Y.z | x > 0) introduce new backward compatible functionality.
  In addition to verifying new features, it's essential to also run integrative tests to confirm that existing functionality is preserved, verifying functional and analytical performance.
- Major versions (X.y.z | X > 0) introduce backward incompatible changes.
  Testing for these releases should thorough, encompassing all levels of testing.

There is specific functionality of nf-core tools that allows developers to easily [update a version of a pipeline when necessary](https://nf-co.re/docs/nf-core-tools/pipelines/bump-version).
Pipelines must be released with stable release tags. Releases must use GitHub releases and keep a detailed changelog file.

Modules used in an nf-core pipeline use fixed software tools inside the container engine (such as docker, singularity, conda). The container packages up the software and all its dependencies so the application runs reliably in any computing environment. In a nf-core module, each of the versions used in a package [must be emitted as per the nf-core guidelines](https://nf-co.re/docs/guidelines/components/modules#emission-of-versions). These versions are documented within a pipeline release.

The [nf-core pipeline releases include checklists to evaluate these aspects](https://nf-co.re/docs/checklists/pipeline_release).

### Code and software development process quality

#### Change Management

In software development, change management refers to the process of tracking and controlling modifications to code and documentation throughout the _software development lifecycle_, ensuring transparency, accountability, and risk mitigation. In regulated bioinformatics environments, change management practices are essential for maintaining compliance, and reproducibility. Within nf-core pipelines, change management is structured as follows:

- Requests for changes, bug reports, and enhancement suggestions can be submitted by any user or community member, ensuring transparent and open community-driven improvement.
- Each pipeline and module follows a Pull Request (PR) template checklist, which helps contributors meet minimum submission requirements.
- Proposed changes must include system and unit tests, which are automatically validated through the continuous integration/continuous deployment (CI/CD) framework, reducing manual testing overhead.
- Changes to development branches require peer review, with each PR needing at least one review before merging into dev branches and at least two reviews before merging into the _main_ branch.
- Automated tests are triggered on each PR to confirm that existing functionality remains unaffected.
- Automated linting checks are performed on each PR, enforcing coding standards and preventing stylistic issues.
- Direct changes to the main branch are not permitted, protecting the integrity of the production-ready code.
- During pipeline release, reviewers must verify that the pipeline adheres to nf-core’s central principles (such as reproducibility, thorough documentation, and compliance with the nf-core template). Any new pipeline submission requires approval from the nf-core core team before integration into the nf-core repository.

**Do I need to re-validate my pipeline every time a change is made?**
In a validated environment, the Risk Assessment process will determine the level of testing required for each change. Minor or non-impactful changes may require testing related to the specific component but major changes may require broader re-validation. Significant changes to the entire pipeline would need a complete re-validation to ensure compliance and integrity.

#### Security

- Patching and updates (including frequency, monitoring of vulnerabilities and third party libraries), requirements management and technical documentation (traceability, reusability, granularity, updates)

### Documentation

#### General documentation

Make sure the pipeline documentation is available and complete.
This should cover at least general aspects of the pipeline and provide a functional default config enabling users to run a basic example.
It should cover subworkflow specific options of the pipeline if there are multiple paths available within a pipeline.

#### Standard operating procedure (SOP)

Establish comprehensive step-by-step instructions, that allow anyone operating the pipeline to do a full run-through for the validated use-case in a consistent way.
It should mention any quality checkpoints or acceptance criteria that need to be applied.

## Testing

#### Functional tests

Functional tests are the tests that nf-core provides to a large extent already for you.
These validate that modules, subworkflows and entire workflow
work functionally, e.g. can be run and produce outputs. These do not cover the full requirements of a validation of an analysis pipeline, which involves integrative tests too.

#### Integrative tests

If you are interested in validating an nf-core pipeline, you are responsible for designing and executing integrative tests that comply with regulatory requirements.
These typically include running validation within your target environment, with data that you will experience during your production setup, e.g. data coming from a sequencing provider using a special kit and in a specific format.
We refer to these tests as integrative tests, which is slightly

TODO integrate this together in one coherent seection

nf-core provides several levels of functional tests for pipelines at each potential stage that composes a pipeline:

- **modules**: We have nf-tests that cover the most atomic units of a pipeline - modules and snapshot the inputs and outputs of a module
- **subworkflows**: We provide nf-tests that cover combined modules (a subworkflow, a certain set of modules within a pipeline)
- **workflows**: We provide test profiles that run the entire pipeline end-to-end with a profile for all available potential subroutes through the analysis pipeline itself

Our plan is to add **analytical performance** tests for pipelines that snapshot and test analytical performance of the pipeline.

We advise you to employ [nf-schema] to perform runtime validation of the config parameters and/or parsed sample sheets.

integrative testing

### Features to include:

Pipeline level:

- Integration testing: an analysis of the test in the production environment with real data
  - Compare the performance of the test system in your dataset with those
    specifications defined by the user. This includes the following performance
    characteristics:
    • Accuracy
    • Precision
    • Reportable range [if applicable]
    • Reference intervals/range (normal values) for the laboratory’s patient population [if applicable]
- Controls to be included to unit tests: [if applicable]
  • Positive control
  • Negative control
  • Additional controls (for example PCR reagent controls, amplification control gene, calibration curve,... )
- Set of expected results for all controls.
- Set assay acceptance criteria
- Set rejection criteria.
  - Add automated quality checks that will be stopping points for the pipeline, if fail.

- Store results to an automated report / stats file
- Automate risk management based on results stored in the stats file

## Maintenance

- Continous development
- Collect bug reports and if possible write a test that covers the affected code.
- Collect testing logs and a history of benchmarking metrics
- Prioritize suggestions for new functionality
- User communication / communication guidelines
- Nf-core template updates create a new minor release at minimum --> not just a patch release
