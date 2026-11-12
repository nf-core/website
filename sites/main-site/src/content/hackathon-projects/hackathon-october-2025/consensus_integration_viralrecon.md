---
title: Viral consensus integration subworkflow using evolutionary priors
category: components
intro_video: ""
slack: https://nfcore.slack.com/archives/C0110R22NH3
leaders:
  GERMAN00VP:
    name: Germán Vallejo
    slack: https://nfcore.slack.com/team/U095H3KMQJX
---

## Project Aim

Current viral consensus generation in nf-core/viralrecon relies on either mapping-based or de novo assembly approaches, each with inherent limitations. Mapping consensuses can miss novel variations in low-coverage regions, while de novo assemblies may produce evolutionarily implausible sequences in highly variable regions.

This project proposes a novel subworkflow that intelligently integrates both mapping and ABACAS-guided de novo consensuses using empirical evolutionary priors. The method employs a sliding-window quality control system that evaluates the evolutionary plausibility of de novo sequences against prior distributions built from large viral sequence databases, ensuring only reliable regions are incorporated into the final consensus.

The implementation is available at [PriorCons repository](https://github.com/GERMAN00VP/PriorCons) and outputs include integrated consensus FASTA files along with comprehensive QC metrics in JSON and TSV formats.

## Goals

### 1. **Subworkflow Implementation**

Create a complete viral consensus integration subworkflow that:

- Performs coverage-based masking of ABACAS consensus with relaxed thresholds (supported by prior information)
- Aligns reference, mapping consensus, and ABACAS consensus sequences
- Optionally builds empirical priors from user-provided databases when priors are not specified
- Executes consensus integration using evolutionary priors for quality control
- Performs final polishing with Pilon using original reads

### 2. **Module Integration**

- Add two new modules to nf-core/modules: `build_priors` and `integrate_consensus`
- Ensure compatibility with existing nf-core/viralrecon workflow structure
- Implement comprehensive testing and documentation

### 3. **Success Metrics**

- **(minimal success)** Functional subworkflow that integrates mapping and de novo consensuses using pre-built priors
- **(success)** Complete implementation including optional prior construction from databases and Pilon polishing
- **(bonus)** Extend the approach to support multiple viral pathogens

## Technical Approach

The subworkflow will leverage the PriorCons methodology which uses:

- **Empirical priors**: Probability distributions built from multiple sequence alignments modeling expected genetic variation
- **Sliding-window QC**: Normalized log-likelihood scoring to evaluate evolutionary plausibility
- **Intelligent merging**: Strategic use of ABACAS sequences to fill missing regions in mapping consensus only when evolutionarily validated

This approach significantly improves consensus quality in highly variable viral genomes while maintaining computational efficiency and providing comprehensive quality metrics.

## Expected Expertise Levels

This project integrates the existing PriorCons tool into nf-core/viralrecon. Work focuses on workflow design, module implementation, and testing.

- **Beginner:** Contribute to documentation, review schema and outputs, and assist with basic workflow tests.
- **Intermediate:** Implement and test the new Nextflow modules (`build_priors`, `integrate_consensus`) and integrate them into the viralrecon structure.
- **Advanced:** Optimize the integration for nf-core standards, ensure compatibility with existing components, and refine parameter settings.

Anyone working with viral sequencing data is welcome to contribute!
