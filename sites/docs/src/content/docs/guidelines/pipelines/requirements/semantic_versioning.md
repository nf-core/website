---
title: Semantic versioning
subtitle: Pipelines must use semantic versioning.
menu:
  main:
    weight: 90
---

- Pipelines must be released with stable release tags.
- Releases must use GitHub releases and include a [detailed changelog](https://keepachangelog.com/en/1.0.0/) file.
- The first release of a pipline should be version `1.0.0` (we strongly discourage pre-releases).
- Release version tags must be numerical only (no `v` prefix) and should follow [semantic versioning](https://semver.org/) rules: `[major].[minor].[patch]`

## Semantic Versioning Criteria

When deciding on the release type, consider the impact on users and whether they need to make changes to their existing workflows or configurations. If in doubt, opt for a higher version increment to clearly communicate the extent of changes to users.

Classify your release according to the highest applicable category in this order: **Major**, **Minor**, then **Patch**. A single match in a higher category overrides all lower categories.

### Major Release (e.g., 1.4.3 to 2.0.0)

Changes that are **not backwards compatible** and may require users to adjust their existing workflows or configurations.

- **Breaking Changes**:
  - Modifications to pipeline user interaction that are not backwards compatible.
  - Renaming or deprecation of existing parameter(s).
  - Changes to input sample sheet specifications (adding, dropping, or renaming mandatory columns).
  - Major alterations to existing output formats (e.g., BAM to CRAM, or changed TSV headers names) or directory structures.
  - Removal of previously supported features or tools.
  - Upgrades to dependencies that introduce breaking changes. For example if the new version of the tool requires additional parameters, removes parameters passed by the pipeline, or changes its outputs in a way that necessitates changes to module output definitions.

:::note
A major release does not necessarily require new functionality; it can result solely from breaking changes.
:::

### Minor Release (e.g., 1.4.3 to 1.5.0)

Introduction of new features that are **backwards compatible** and do not alter existing functionality.

- **Feature Additions**:
  - Adding new parameters that introduce additional functionality or options to existing modules (without fundamentally changing output file types or structure).
  - Incorporating new optional or default processes without changing user interaction with previous steps.
  - Introducing new output files or metrics without modifying existing ones.
  - Supporting new input formats or data types without removing support for old ones.
- **Improvements**:
  - Implementing performance enhancements that don't affect the pipeline's interface (e.g. refined default memory or CPU requirements).
  - Expanding documentation or adding new examples without altering core functionality.
  - A template update to a newer nf-core template.

### Patch Release (e.g., 1.4.3 to 1.4.4)

Minor changes that do not affect user interaction or existing functionality.

- **Bug Fixes**:
  - Correcting minor bugs that do not change how the user interacts with the pipeline.
  - Fixing the logic of pipeline operators to ensure correct inputs for a process.
  - Inserting missing arguments into process command lines that should have existed before (but without modifying output file types/output file content structure).
- **Updates and Corrections**:
  - Updating a container build version of a tool without changing the tool's version itself (e.g. so it can execute correctly).
  - Correcting typos in documentation or code comments.
  - Updating metadata or non-functional elements (e.g., README, CONTRIBUTING files).
  - Making minor adjustments to log messages or error handling.
  - Applying small optimizations that don't affect functionality or user interaction.
