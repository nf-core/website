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

**Versioning Examples**

Starting from release version `1.4.3`, bumping the version to:

## Release type checklist

If your release matches any point in the checklists below, classify it according to the highest applicable category in this order: major, minor, then patch. A single match in a higher category overrides all lower categories.

### Major Release (e.g., 1.0.0 to 2.0.0)

- Changes that are not backwards compatible in terms of pipeline user interaction
- Renaming or deprecation of existing parameter(s)
- Changes to input sample sheet specification (adding, dropping, renaming of mandatory columns)
- Significant changes to the pipeline's core functionality or workflow
- Major changes to output formats or directory structures
- Removal of previously supported features or tools
- Changes that require users to modify their existing configuration files or scripts
- Upgrades to dependencies that introduce breaking changes

## Minor Release (e.g., 1.4.0 to 1.5.0)

- Addition of new features without changing existing functionality
- Adding a new parameter that introduces new functionality/options to an existing module
- Adding a new optional or default process that doesn't change user interaction with previous steps
- Introducing new output files or metrics without altering existing ones
- Adding support for new input formats or data types
- Implementing performance improvements that don't affect the pipeline's interface
- Expanding documentation or adding new examples without changing core functionality
- Adding new optional dependencies or tools

## Patch Release (e.g., 1.4.3 to 1.4.4)

- Minor fixes such as bug corrections that do not modify how the user interacts with the pipeline
- Inserting a missing argument to a process command line
- Updating a container build version of a tool that is otherwise the same version of the tool itself
- Fixing the logic of pipeline operators to ensure correct inputs for a process
- Correcting typos in documentation or code comments
- Updating metadata or non-functional elements (e.g., README, CONTRIBUTING files)
- Making minor adjustments to log messages or error handling
- Applying small optimizations that don't affect functionality or user interaction

When deciding on the release type, it should. be considered the impact on users and whether they need to make changes to their existing workflows or configurations. If in doubt, it's generally safer to opt for a higher version increment to clearly communicate the extent of changes to users.
  - Fixing the logic of the pipeline operators to ensure correct inputs for a process
- **`1.5.0`**: A minor release that adds new features without changing existing functionality. Examples:
  - Adding a new parameter that adds new functionality/options to an existing module
  - Add a new optional or default process that does not change the way the user interacts with previous steps of a pipeline
- **`2.0.0`**: A major release where execution interaction, inputs or result structures are no longer backwards compatible in terms of pipeline user interaction. Examples:
  - Renaming or deprecations of existing parameter(s)
  - Change to input sample sheet specification (adding, dropping, renaming of mandatory columns) or change the input samplesheet specification
  - Note that a major release DOES NOT necessarily require new functionality
