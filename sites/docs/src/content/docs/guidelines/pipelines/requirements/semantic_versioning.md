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

- **`1.4.4`**: A patch release for minor fixes such as bug corrections that do not modify how the user interacts with the pipeline.
- **`1.5`**: A minor release that adds new features without changing existing functionality.
- **`2.0`**: A major release where inputs or results are no longer backwards compatible, or when you rename or deprecate a parameter, or change the input samplesheet specification. This does not require major changes in workflow functionality—alterations to the pipeline's 'interface' are sufficient to justify a new major version.
