---
title: Overview
subtitle: Overview on regulatory aspects for nf-core pipelines
weight: 1
---

# Introduction

nf-core aims to produce high-quality Nextflow pipelines that facilitate validation and regulatory compliance through a comprehensive metrics-based approach. The community is organized around the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory), where stakeholders collaborate to define and implement validation standards.

## Document Version

2.0.0 draft

## nf-core Validation Strategy

nf-core has developed a systematic approach to pipeline validation that balances community-driven quality assurance with institutional flexibility:

### Our Approach

- **Metrics Collection**: nf-core automatically collects and maintains comprehensive quality metrics for all pipelines through the nf-core/stats service
- **Validation Readiness Reports**: Each pipeline release includes an automated PDF report containing validation-relevant metrics and quality indicators
- **Community Standards**: The nf-core community maintains consistent quality standards across all pipelines through established guidelines and testing frameworks
- **User Flexibility**: Institutions use nf-core reports as a foundation for their own validation processes, tailored to specific regulatory requirements and infrastructure

### What nf-core Provides

✅ **Comprehensive metrics collection** via nf-core/stats
✅ **Automated validation readiness reports** with each pipeline release
✅ **Quality assurance frameworks** (testing, documentation, versioning)
✅ **Implementation guidance** for institutional validation
✅ **Continuous improvement** through gap identification and community feedback

### What Users Are Responsible For

🔄 **Institution-specific validation** within their target environment
🔄 **Regulatory compliance** according to applicable standards
🔄 **Infrastructure validation** (CSV and system-level requirements)
🔄 **Risk assessment** for their specific use cases
🔄 **Ongoing maintenance** and re-validation as needed

> [!NOTE]
> nf-core provides the foundation for validation through standardized metrics and quality assurance, but cannot provide complete "off-the-shelf" validation reports. Each institution must perform their own validation according to their specific regulatory requirements and risk assessment.

## nf-core Metrics & Reporting System

### nf-core/stats Service

The [nf-core/stats](https://github.com/nf-core/stats) service serves as the central repository for collecting, storing, and analyzing pipeline quality metrics. This system automatically:

- **Tracks pipeline metrics** across all nf-core pipelines and releases
- **Monitors community engagement** (contributors, users, issue resolution)
- **Assesses code quality** (testing coverage, documentation completeness)
- **Evaluates pipeline maturity** (release history, stability indicators)
- **Provides historical data** for trend analysis and continuous improvement

### Automated Validation Reports

Each nf-core pipeline release automatically generates a comprehensive PDF validation readiness report containing:

#### Pipeline Overview
- Release information and semantic versioning details
- Community metrics (contributors, users, maintenance activity)
- Licensing and governance information

#### Quality Metrics
- **Testing Coverage**: Functional and integration test results
- **Documentation Quality**: Completeness and clarity assessments
- **Code Quality**: Linting, standards compliance, security scanning
- **Dependency Management**: Tool versions, container specifications

#### Community Indicators
- **Maintenance Activity**: Issue resolution rates, update frequency
- **User Adoption**: Download statistics, citation metrics
- **Contributor Engagement**: Review processes, contribution patterns

#### Risk Assessment Factors
- **Complexity Analysis**: Pipeline complexity scoring
- **Stability Indicators**: Error rates, failure patterns
- **Change Management**: Update frequency and impact analysis

### Gap Identification and Improvement

The nf-core community continuously evaluates and enhances the metrics collection system:

- **Community Feedback**: Regular review of validation requirements across different domains
- **Standards Evolution**: Adaptation to changing regulatory landscapes
- **Metric Enhancement**: Addition of new quality indicators based on community needs
- **Tooling Development**: Improvement of automated collection and reporting systems

> [!TIP]
> To request additional metrics or report enhancements, please engage with the [regulatory special interest group](https://nf-co.re/special-interest-groups/regulatory) or open an issue in the [nf-core/stats repository](https://github.com/nf-core/stats).

## Regulatory Context

### Risk-Based Validation Principles

Pipeline validation follows risk-based approaches that consider the development, implementation, and integration of analysis pipelines as potential sources of risk. Common risks include:

- **Execution Failures**: Risk of pipeline not producing expected outputs
- **Data Integrity**: Risk of data corruption or loss during processing
- **Analytical Accuracy**: Risk of incorrect or inconsistent results
- **Reproducibility**: Risk of inability to replicate results across runs
- **Security**: Risk of unauthorized access or data breaches

These risks are mitigated through appropriate measures such as comprehensive testing, documentation, version control, and quality assurance processes.

### Applicable Standards

Regulatory requirements vary by jurisdiction and intended use. Common standards that may apply to bioinformatics pipelines include:

#### Software as Medical Device (SaMD)
- [FDA SaMD Guidance](https://www.fda.gov/medical-devices/digital-health-center-excellence/software-medical-device-samd)
- [SaMD General Guidelines](https://www.greenlight.guru/blog/samd-software-as-a-medical-device)

#### Laboratory Developed Tests (LDT)
- [FDA LDT Regulations](https://www.fda.gov/medical-devices/in-vitro-diagnostics/laboratory-developed-tests)
- [CLIA Validation Requirements](https://www.cms.gov/regulations-and-guidance/legislation/clia/downloads/6064bk.pdf)

#### Medical Device Registration
- [CE Mark Registration](https://europa.eu/youreurope/business/product-requirements/labels-markings/ce-marking/index_en.htm)
- [EU Medical Device Guidelines](https://health.ec.europa.eu/system/files/2021-10/mdcg_2021-24_en_0.pdf)

#### Harmonization Efforts
- [FDA Consensus Standards](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/appropriate-use-voluntary-consensus-standards-premarket-submissions-medical-devices)
- [EU Harmonized Standards](https://single-market-economy.ec.europa.eu/single-market/european-standards/harmonised-standards_en)

> [!NOTE]
> Computerized Systems Validation (CSV) requirements for infrastructure and IT systems are outside the scope of this guidance. Please consult with your IT and compliance teams for CSV requirements.

## Using nf-core Validation Reports

### Interpreting the Automated Reports

The validation readiness reports provide a standardized foundation for institutional validation. Key sections to review include:

#### 1. **Executive Summary**
- Overall pipeline maturity score
- Critical risk indicators
- Recommended validation approach

#### 2. **Quality Metrics Analysis**
- Testing coverage and pass rates
- Documentation completeness scores
- Code quality indicators
- Community engagement metrics

#### 3. **Risk Assessment Matrix**
- Identified risk factors specific to the pipeline
- Severity and likelihood assessments
- Recommended mitigation strategies

#### 4. **Compliance Readiness**
- Alignment with common regulatory requirements
- Gap analysis for specific use cases
- Documentation and traceability indicators

### Institution-Specific Validation Process

Using the nf-core report as a foundation, institutions should follow these steps:

#### Phase 1: Initial Assessment
1. **Review nf-core validation report** for the target pipeline version
2. **Define intended use** and applicable regulatory requirements
3. **Conduct risk assessment** based on institutional context
4. **Identify validation scope** and required evidence

#### Phase 2: Gap Analysis
1. **Compare nf-core metrics** against institutional requirements
2. **Identify additional testing** needs for specific use cases
3. **Assess infrastructure requirements** and CSV considerations
4. **Plan integration testing** with institutional data and systems

#### Phase 3: Validation Execution
1. **Perform functional verification** in target environment
2. **Execute integration testing** with representative datasets
3. **Validate analytical performance** against established criteria
4. **Document evidence** of compliance and performance

#### Phase 4: Ongoing Maintenance
1. **Monitor pipeline updates** and assess impact
2. **Maintain validation status** through change control
3. **Update validation** as requirements evolve
4. **Participate in community** feedback and improvement

### Infrastructure Considerations

#### Computing Environment
- **Reproducibility**: Ensure consistent execution across runs
- **Security**: Implement appropriate access controls and audit trails
- **Performance**: Validate pipeline performance under expected loads
- **Backup/Recovery**: Establish data protection and recovery procedures

#### Data Management
- **Input Validation**: Verify data integrity and format compliance
- **Output Verification**: Confirm expected outputs and quality metrics
- **Traceability**: Maintain complete audit trails of data processing
- **Retention**: Implement appropriate data retention policies

> [!IMPORTANT]
> The nf-core validation report provides a foundation, but each institution must validate the pipeline within their specific environment, with their data, and according to their regulatory requirements.

## Community Standards and Guidelines

### Pipeline Quality Requirements

All nf-core pipelines must meet standardized quality criteria that support validation efforts:

#### Development Standards
- **Semantic Versioning**: Clear version management with defined backward compatibility
- **Code Quality**: Adherence to coding standards and automated linting
- **Testing Framework**: Comprehensive functional and integration tests
- **Documentation**: Complete usage instructions and parameter documentation

#### Change Management
- **Pull Request Process**: Mandatory peer review for all changes
- **Automated Testing**: CI/CD pipeline validation for every change
- **Release Process**: Structured release workflow with quality checkpoints
- **Issue Tracking**: Transparent bug reporting and feature request management

#### Security and Compliance
- **Dependency Management**: Fixed software versions in containerized environments
- **Vulnerability Scanning**: Regular security assessment of dependencies
- **License Compliance**: Clear licensing for all components
- **Access Control**: Controlled access to pipeline repositories and releases

### Testing Framework

nf-core provides comprehensive testing at multiple levels to ensure pipeline quality:

#### Functional Testing
- **Module Tests**: Validate individual pipeline components using [nf-test](https://www.nf-test.com/)
- **Subworkflow Tests**: Test combined modules and workflow sections
- **Pipeline Tests**: End-to-end testing with multiple test profiles
- **Regression Testing**: Automated testing to prevent functional regressions

#### Performance Testing
- **Resource Usage**: Monitor computational requirements and efficiency
- **Scalability**: Test pipeline performance across different dataset sizes
- **Benchmarking**: Track performance metrics over time

#### Integration Testing
- **Environment Testing**: Validate execution across different computing environments
- **Data Format Testing**: Ensure compatibility with various input formats
- **Output Validation**: Verify expected outputs and quality metrics

### Documentation Standards

Comprehensive documentation supports validation and ensures consistent usage:

#### User Documentation
- **Usage Instructions**: Clear step-by-step execution guidance
- **Parameter Documentation**: Complete description of all configurable options
- **Output Documentation**: Description of all pipeline outputs and their interpretation
- **Troubleshooting Guides**: Common issues and resolution strategies

#### Technical Documentation
- **Architecture Overview**: Pipeline structure and workflow design
- **Module Documentation**: Individual component specifications
- **Dependency Documentation**: Complete list of software dependencies and versions
- **Configuration Guidelines**: Environment setup and configuration requirements

## Implementation Guidance

### Step-by-Step Validation Process

#### 1. Preparation Phase
```
□ Download pipeline validation report from nf-core release
□ Define intended use and regulatory requirements
□ Identify applicable standards and compliance frameworks
□ Assemble validation team with appropriate expertise
□ Establish validation protocols and acceptance criteria
```

#### 2. Risk Assessment
```
□ Review nf-core risk indicators and community metrics
□ Conduct institution-specific risk analysis
□ Identify critical control points and testing requirements
□ Document risk mitigation strategies
□ Establish ongoing monitoring procedures
```

#### 3. Functional Validation
```
□ Set up validated computing environment
□ Execute nf-core test profiles successfully
□ Perform integration testing with institutional data
□ Validate analytical performance against specifications
□ Document all testing evidence and results
```

#### 4. Operational Validation
```
□ Establish standard operating procedures (SOPs)
□ Train personnel on validated procedures
□ Implement quality control measures
□ Establish data management and retention procedures
□ Create incident response and change control processes
```

#### 5. Ongoing Maintenance
```
□ Monitor pipeline updates and security advisories
□ Assess impact of changes through change control
□ Maintain validation status through periodic review
□ Update procedures as requirements evolve
□ Participate in community feedback and improvement
```

### Gap Analysis and Improvement

#### Identifying Validation Gaps
- **Regulatory Requirements**: Compare nf-core metrics against specific compliance needs
- **Technical Requirements**: Assess infrastructure and integration requirements
- **Analytical Requirements**: Evaluate performance specifications and acceptance criteria
- **Documentation Requirements**: Identify additional documentation needs

#### Improvement Pathways
- **Community Engagement**: Participate in nf-core regulatory working groups
- **Metric Enhancement**: Propose additional metrics through nf-core/stats
- **Documentation Contributions**: Contribute to pipeline documentation and guides
- **Testing Contributions**: Develop and share additional test cases

### Maintenance and Continuous Improvement

#### Change Management
- **Version Updates**: Systematic evaluation of pipeline updates
- **Impact Assessment**: Risk-based analysis of changes
- **Re-validation**: Appropriate level of testing for changes
- **Documentation Updates**: Maintain current validation documentation

#### Quality Monitoring
- **Performance Tracking**: Monitor pipeline performance over time
- **Issue Resolution**: Track and resolve validation-related issues
- **Trend Analysis**: Analyze patterns in pipeline usage and performance
- **Continuous Improvement**: Regular review and enhancement of validation processes

> [!SUCCESS]
> By following this structured approach and leveraging nf-core's comprehensive metrics and reporting system, institutions can efficiently validate nf-core pipelines while maintaining the flexibility to meet their specific regulatory and operational requirements.

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
