---
title: Semantic versioning
subtitle: Pipelines must use semantic versioning.
menu:
  main:
    weight: 90
---

- Pipelines must be released with stable release tags.
- Releases must use GitHub releases and include a [detailed changelog](https://keepachangelog.com/en/1.0.0/) file.
- Release version tags must be numerical only (no `v` prefix) and should follow [semantic versioning](https://semver.org/) rules: `[major].[minor].[patch]`

**Versioning Examples**

Starting from release version `1.4.3`, bumping the version to:

- **`1.4.4`**: A patch release for minor fixes such as bug corrections that do not modify how the user interacts with the pipeline. Examples:
  - Inserting a missing argument to a process command line
  - Updating a container build version of a tool that is otherwise the same version of the tool itself
- **`1.5`**: A minor release that adds new features without changing existing functionality. Examples:
  - Adding a new parameter that adds new functionality/options to an existing module
  - Add a new optional or default process that does not change the way the user interacts with previous steps of a pipeline
- **`2.0`**: A major release where execution interaction, inputs or results are no longer backwards in terms of pipeline user interaction. Examples:
  - Renaming or deprecations of existing parameter(s)
  - Change to input sample sheet specification (adding, dropping, renaming of mandatory columns) or change the input samplesheet specification.
  - Note that a major release DOES NOT necessarily require new functionality
