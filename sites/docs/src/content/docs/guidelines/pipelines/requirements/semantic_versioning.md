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
  - Major alterations to output formats or directory structures.
  - Removal of previously supported features or tools.
  - Upgrades to dependencies that introduce breaking changes.
- **Significant Functional Changes**:
  - Overhauls to the pipeline's core functionality or workflow.

**Note:** A major release does not necessarily require new functionality; it can result solely from breaking changes.

### Minor Release (e.g., 1.4.3 to 1.5.0)

Introduction of new features that are **backwards compatible** and do not alter existing functionality.

- **Feature Additions**:
  - Adding new parameters that introduce additional functionality or options to existing modules.
  - Incorporating new optional or default processes without changing user interaction with previous steps.
  - Introducing new output files or metrics without modifying existing ones.
  - Supporting new input formats or data types.
- **Improvements**:
  - Implementing performance enhancements that don't affect the pipeline's interface.
  - Expanding documentation or adding new examples without altering core functionality.
  - Adding new optional dependencies or tools.

### Patch Release (e.g., 1.4.3 to 1.4.4)

Minor changes that do not affect user interaction or existing functionality.

- **Bug Fixes**:
  - Correcting minor bugs that do not change how the user interacts with the pipeline.
  - Fixing the logic of pipeline operators to ensure correct inputs for a process.
  - Inserting missing arguments into process command lines.
- **Updates and Corrections**:
  - Updating a container build version of a tool without changing the tool's version itself.
  - Correcting typos in documentation or code comments.
  - Updating metadata or non-functional elements (e.g., README, CONTRIBUTING files).
  - Making minor adjustments to log messages or error handling.
  - Applying small optimizations that don't affect functionality or user interaction.
