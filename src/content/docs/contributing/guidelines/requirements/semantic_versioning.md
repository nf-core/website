---
title: Semantic versioning
subtitle: Pipelines must use semantic versioning.
menu:
  main:
    weight: 90
---

Pipelines must be released with stable release tags.
Releases must use GitHub releases and [keep a detailed changelog](https://keepachangelog.com/en/1.0.0/) file.

Release version tags must be numerical only (no `v` prefix) and should follow [semantic versioning](https://semver.org/) rules: `[major].[minor].[patch]`

For example, starting with with a release version `1.4.3`, bumping the version to:

- `1.4.4` would be a patch release for minor things such as fixing bugs.
- `1.5` would be a minor release, for example adding some new features.
- `2.0` would correspond to the _major_ release where results would no longer be backwards compatible.
