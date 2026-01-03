---
title: Pipelines
subtitle: Specifications for nf-core pipelines
weight: 1
parentWeight: 20
---

nf-core pipelines are end-to-end analysis workflows that combine multiple tools and processes to analyze biological data. These specifications define the standards for developing high-quality, reproducible pipelines that can be shared and maintained by the nf-core community.

The specifications are organized into **requirements** (mandatory standards all nf-core pipelines MUST follow) and **recommendations** (best practices pipelines SHOULD adopt where applicable). While these specifications are mandatory for pipelines contributed to the nf-core repository, they represent proven best practices for developing professional-grade Nextflow pipelines that can benefit any workflow development project.

:::tip
If you're thinking of adding a new pipeline to nf-core or contributing to an existing pipeline, see [adding a new pipeline](/docs/tutorials/adding_a_pipeline/overview) and [contributing to a pipeline](/docs/tutorials/contributing_to_nf-core/contributing_to_pipelines).
:::

<!-- TODO: Update links -->

:::warning{title="Build with the community"}
nf-core is a community of builders, not a registry of pipelines. Propose your pipeline and develop with us, rather than proposing an already-complete pipeline. See recommendation: [Build with the community](/docs/guidelines/pipelines/recommendations/build_with_community)
:::

## Requirements

The following specifications define mandatory standards that all nf-core pipelines MUST follow:

- **[Nextflow](/docs/guidelines/pipelines/requirements/nextflow):** Pipelines MUST be built using Nextflow and adhere to DSL2 syntax standards.
- **[Community owned](/docs/guidelines/pipelines/requirements/community_owned):** Pipelines are owned collectively by the community with contributions welcomed from all members.
- **[Identity and branding](/docs/guidelines/pipelines/requirements/identity_branding):** Primary development MUST occur within the nf-core GitHub organisation with appropriate branding.
- **[Workflow specificity](/docs/guidelines/pipelines/requirements/workflow_specificity):** Each data or analysis type MUST have only one canonical nf-core pipeline to avoid fragmentation.
- **[Workflow size](/docs/guidelines/pipelines/requirements/workflow_size):** Pipelines MUST have appropriate scope—comprehensive enough to be useful but focused enough to maintain.
- **[Workflow name](/docs/guidelines/pipelines/requirements/workflow_name):** Pipeline names MUST be lowercase, use hyphens for separation, and avoid punctuation.
- **[Use the template](/docs/guidelines/pipelines/requirements/use_the_template):** All pipelines MUST be built using the nf-core template and maintain synchronization through automated processes.
- **[Software license](/docs/guidelines/pipelines/requirements/mit_license):** Pipelines MUST be open source and released under the MIT license.
- **[Bundled documentation](/docs/guidelines/pipelines/requirements/docs):** Comprehensive documentation MUST be bundled with the pipeline and hosted on the nf-core website.
- **[Docker support](/docs/guidelines/pipelines/requirements/docker):** Software dependencies MUST be containerized using Docker with specific version tags.
- **[Continuous integration testing](/docs/guidelines/pipelines/requirements/ci_testing):** Pipelines MUST implement automated CI tests to verify functionality.
- **[Semantic versioning](/docs/guidelines/pipelines/requirements/semantic_versioning):** Pipelines MUST follow semantic versioning with stable release tags.
- **[Standardised parameters](/docs/guidelines/pipelines/requirements/parameters):** Pipelines MUST use standardized parameter names and conventions for consistency across nf-core.
- **[Single command](/docs/guidelines/pipelines/requirements/single_command):** Pipelines MUST be executable with a single command without requiring manual intervention.
- **[Keywords](/docs/guidelines/pipelines/requirements/keywords):** Pipelines MUST have comprehensive documentation and appropriate GitHub repository keywords for discoverability.
- **[Pass lint tests](/docs/guidelines/pipelines/requirements/linting):** Pipelines MUST pass all `nf-core pipelines lint` tests without failures.
- **[Credits and Acknowledgements](/docs/guidelines/pipelines/requirements/acknowledgements):** Pipelines MUST properly acknowledge and credit prior work, tools, and contributors.
- **[Minimum inputs](/docs/guidelines/pipelines/requirements/minimum_inputs):** Pipelines MUST be designed to run with minimal required input, using sensible defaults.
- **[Use nf-core git branches](/docs/guidelines/pipelines/requirements/git_branches):** Pipelines MUST use standardized branch structure including `master`, `dev`, and `TEMPLATE`.

## Recommendations

The following specifications define best practices that all nf-core pipelines SHOULD follow where possible or appropriate:

- **[Use Bioconda](/docs/guidelines/pipelines/recommendations/bioconda):** Software dependencies SHOULD be packaged using Bioconda and distributed through BioContainers for reproducibility.
- **[File formats](/docs/guidelines/pipelines/recommendations/file_formats):** Pipelines SHOULD use community-accepted modern file formats such as CRAM for better compression and compatibility.
- **[Testing](/docs/guidelines/pipelines/recommendations/testing):** Pipelines SHOULD use nf-test to verify successful completion with valid outputs using minimal test datasets.
- **[DOIs](/docs/guidelines/pipelines/recommendations/dois):** Pipeline releases SHOULD be assigned digital object identifiers (DOIs) through Zenodo for citability.
- **[Cloud compatible](/docs/guidelines/pipelines/recommendations/cloud_compatible):** Pipelines SHOULD be tested on cloud computing environments to ensure portability.
- **[Publication credit](/docs/guidelines/pipelines/recommendations/publication_credit):** Pipeline publications SHOULD acknowledge the nf-core community and cite the nf-core framework paper.
- **[Build with the community](/docs/guidelines/pipelines/recommendations/build_with_community):** Pipeline development SHOULD occur collaboratively with the community from the beginning rather than proposing complete pipelines.
- **[Custom Docker images](/docs/guidelines/pipelines/recommendations/custom_containers):** Custom containers not available through BioContainers SHOULD be mirrored on quay.io for long-term reproducibility.
