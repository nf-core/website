---
title: Create a ROCrate
subtitle: Create a Research Object (RO) Crate
shortTitle: rocrate
---

The `nf-core pipelines rocrate` command generates a Research Object (RO) Crate for your pipeline.
RO Crates are an open-source standard which in this case describe a pipeline and its components in a structured way and help with automated provenance tracking.

You can generate and update a RO Crate using the `nf-core pipelines rocrate` command (which is automatically run with `nf-core pipelines create` and `nf-core pipelines sync`).

This command takes a pipeline directory and attempts to generate a RO Crate for it, relying heavily on the [repo2rocrate library](https://github.com/crs4/repo2rocrate).

:::note
The `author` field is populated based on the git contributors.
The ORCID ID is currently identified by searching the ORCID API for a unique match of the contributor's full name provided from GitHub.
:::
