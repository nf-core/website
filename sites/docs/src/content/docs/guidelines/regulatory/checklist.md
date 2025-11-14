---
title: Pipeline validation checklist
subtitle: A practical checklist for validating nf-core pipelines in your institution
weight: 10
markdownPlugin: checklist
---

This checklist guides you through the validation process for nf-core pipelines based on the principles outlined in the overview documentation. Complete each section according to your institutional requirements and regulatory context.

:::note
This checklist serves as a comprehensive framework, though certain sections vary in detail by design. Some aspects, such as Computer Systems Validation (CSV), are deliberately kept concise as they fall outside the scope of the nf-core community and require institution-specific implementation tailored to individual organizational requirements and regulatory contexts.
:::

# Phase 1: Initial assessment

## Preparation

- [ ] Download the nf-core validation readiness report for your target pipeline version
- [ ] Define the intended use of the pipeline in your institution
- [ ] Identify applicable regulatory requirements (e.g., FDA SaMD, CLIA, CE Mark, LDT)
- [ ] Assemble validation team with appropriate expertise
- [ ] Establish validation protocols and acceptance criteria

## Risk assessment

- [ ] Review nf-core risk indicators from the validation report
- [ ] Conduct institution-specific risk analysis for your use case
- [ ] Identify critical control points (execution failures, data integrity, analytical accuracy, reproducibility, security)
- [ ] Document risk mitigation strategies
- [ ] Define ongoing monitoring procedures

# Phase 2: Gap analysis

## Review nf-core metrics

- [ ] Examine testing coverage and pass rates
- [ ] Review documentation completeness scores
- [ ] Check code quality indicators
- [ ] Assess community engagement metrics (contributors, maintenance activity, issue resolution)

## Identify additional requirements

- [ ] Compare nf-core metrics against institutional requirements
- [ ] Identify additional testing needs for your specific use case
- [ ] Assess infrastructure requirements and CSV considerations
- [ ] Plan integration testing with institutional data and systems
- [ ] Document validation scope and required evidence

# Phase 3: Validation execution

## Computing environment

- [ ] Set up validated computing environment
- [ ] Ensure consistent execution across runs (reproducibility)
- [ ] Implement appropriate access controls and audit trails (security)
- [ ] Validate pipeline performance under expected loads
- [ ] Establish data protection and recovery procedures

## Functional validation

- [ ] Execute all nf-core test profiles successfully
- [ ] Verify module-level tests pass using nf-test
- [ ] Verify subworkflow-level tests pass
- [ ] Validate end-to-end pipeline execution with test profiles
- [ ] Document all testing evidence and results

## Integration testing

- [ ] Execute integration testing with institutional data in production environment
- [ ] Test with data in expected formats from your sequencing provider
- [ ] Validate analytical performance against your specifications (accuracy, precision, reportable range)
- [ ] Include appropriate controls (positive control, negative control, additional controls if applicable)
- [ ] Define and verify expected results for all controls

## Data management

- [ ] Validate input data integrity and format compliance
- [ ] Verify expected outputs and quality metrics
- [ ] Establish complete audit trails for data processing
- [ ] Implement data retention policies
- [ ] Use nf-schema for runtime validation of config parameters and sample sheets
- [ ] Ensure pipelines are running with the `nf-prov` plugin enabled, creating further traceability outputs for auditing purposes

## Quality control

- [ ] Set assay acceptance criteria
- [ ] Set rejection criteria
- [ ] If applicable, add automated quality checks as pipeline stopping points if criteria fail
- [ ] Store results to automated report/stats file
- [ ] Implement automated risk management based on stored results

# Phase 4: Operational validation

## Standard operating procedures

- [ ] Create comprehensive step-by-step instructions for consistent pipeline operation
- [ ] Document validated use-case run-through procedures
- [ ] Define quality checkpoints and acceptance criteria
- [ ] Establish incident response procedures
- [ ] Create change control processes

## Personnel and training

- [ ] Train personnel on validated procedures
- [ ] Establish competency assessment protocols
- [ ] Document training records

## Documentation

- [ ] Verify pipeline documentation is complete
- [ ] Ensure functional default config enables basic example execution
- [ ] Document subworkflow-specific options if multiple paths exist
- [ ] Maintain validation documentation package

# Phase 5: Ongoing maintenance

## Change management

- [ ] Monitor pipeline updates and security advisories from nf-core
- [ ] Assess impact of changes through institutional change control process
- [ ] Perform risk-based analysis for each pipeline update
- [ ] Determine re-validation requirements (component-level, broader, or complete)
- [ ] Track and resolve validation-related issues
- [ ] Update validation documentation as needed

## Quality monitoring

- [ ] Collect and review testing logs
- [ ] Maintain history of benchmarking metrics
- [ ] Collect and prioritize bug reports
- [ ] Write tests covering reported bugs
- [ ] Monitor pipeline performance over time
- [ ] Analyze trends in pipeline usage and performance

## Community engagement

- [ ] Participate in nf-core regulatory special interest group
- [ ] Provide feedback on validation reports and metrics
- [ ] Propose additional metrics through nf-core/stats if needed
- [ ] Contribute to pipeline documentation and test cases
- [ ] Stay informed on nf-core template updates (trigger new minor release minimum)

# Additional considerations

## Versioning and releases

- [ ] Verify pipeline uses semantic versioning
- [ ] Track pipeline version used in validation
- [ ] Monitor backward compatibility across versions
- [ ] Document version-specific validation status

## Licensing and governance

- [ ] Review pipeline license compliance (MIT license)
- [ ] Verify dependency licenses
- [ ] Understand nf-core governance model
- [ ] Review access control to pipeline repositories

## Security

- [ ] Monitor vulnerabilities in dependencies
- [ ] Review container security scans (Seqera Containers)
- [ ] Plan patching and update schedule
- [ ] Assess third-party library risks

:::warning{title="Caution"}
This checklist provides a framework based on nf-core's validation approach. Adapt it to your specific regulatory requirements, institutional policies, and risk assessment outcomes. This is not necessarily complete or applies to all potential cases, but should serve as a go-to-resource that facilitates doing a pipeline validation at your respective institution and also facilitate finding appropriate lines of thought by regulatory authorities for example.
:::
