---
title: List pipelines
subtitle: Listing nf-core pipelines, local and remote
shortTitle: list
weight: 10
parentWeight: 10
---

The `nf-core pipelines list` command shows all available nf-core pipelines with their latest version, publication date, and when you last pulled the pipeline locally.

Example output:

<!-- RICH-CODEX head: 19 -->

![`nf-core pipelines list`](../../../../../assets/images/tools/nf-core-list.svg)

Supply one or more keywords to filter pipelines by titles, descriptions, and topics:

![`nf-core pipelines list rna rna-seq`](../../../../../assets/images/tools/nf-core-list-rna.svg)

Sort results by latest release (`-s release`, default), when you last pulled a local copy (`-s pulled`), alphabetically (`-s name`), or by number of GitHub stars (`-s stars`).

<!-- RICH-CODEX head: 18 -->

![`nf-core pipelines list -s stars`](../../../../../assets/images/tools/nf-core-list-stars.svg)

To return results as JSON output for downstream use, use the `--json` flag.

Archived pipelines are not returned by default. To include them, use the `--show_archived` flag.
