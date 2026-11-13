---
title: Use the template
subtitle: All nf-core pipelines must be built using the nf-core template.
menu:
  main:
    weight: 40
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

All nf-core pipelines MUST be built using the nf-core template.

:::note
If you plan to contribute your pipeline to nf-core, propose it to the community before you create it.
See [Propose before you build](../../../contributing/contribute-new-pipelines#propose-before-you-build).
:::

Workflows MUST be started using the `nf-core pipelines create` command which makes a new git repository and the initial commits and branches.
This ensures that the automated sync process that keeps all nf-core pipelines up to date can work.
See the [sync docs](/docs/tutorials/sync/overview) for details.

:::info{title="Rationale" collapse}
The automated sync process relies on specific git history and branch structures created by the `nf-core pipelines create` command.
Without this initial setup, the sync tool cannot properly track and merge template updates.
:::

Workflow authors SHOULD follow nf-core conventions for filenames and code locations where possible.
