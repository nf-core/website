---
title: Create a ROCrate
subtitle: Create a Research Object (RO) Crate
shortTitle: rocrate
weight: 60
---

The `nf-core pipelines rocrate` command generates a Research Object (RO) Crate for your pipeline.
RO Crates are an open-source standard that describe a pipeline and its components in a structured way and help with automated provenance tracking.

Generate and update a RO Crate using the `nf-core pipelines rocrate` command. This runs automatically with `nf-core pipelines create` and `nf-core pipelines sync`.

The command takes a pipeline directory and generates a RO Crate for it using the [repo2rocrate library](https://github.com/crs4/repo2rocrate).

:::note
The `author` field is populated based on the git contributors.
The ORCID ID is currently identified by searching the ORCID API for a unique match of the contributor's full name provided from GitHub.
:::
