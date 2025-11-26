---
title: Review checklist for nf-core/tools
subtitle: A checklist for reviewing PRs in the nf-core/tools repository
shortTitle: nf-core/tools PR
markdownPlugin: checklist
---

<!-- TODO: Add links to other pages and guide where possible -->

When you review a PR, you examine changes proposed to the `nf-core/tools` repository. You provide constructive feedback on those changes before the maintainers merge them. Your review ensures that the code meets the coding standards of the project, maintains consistency and achieves high quality.

While the infrastructure team oversees the PR review process for nf-core/tools, these guidelines help you review PRs consistently and effectively as a community member.

## General checklist items

- [ ] Check that tests cover the added code (either old or newly added tests).
- [ ] Verify documentation exists or the PR adds it, including checking for corresponding PR updates to nf-core/website if needed.
- [ ] Ensure the code does not hardcode "nf-core"; check it uses `organisation` variable instead.
- [ ] Check the code does not hardcode `modules`/`subworkflows`; verify it uses `component_type` instead.

## Python styling checklist

- [ ] Prefer `pathlib` over `os` module when possible for file operations.
- [ ] Check that the code uses typing.
- [ ] Verify comments are sufficient and clear for readability.
- [ ] Update docstrings as they generate API documentation automatically.
- [ ] If the PR adds a new command, check that the API documentation skeleton includes it.

## Template modifications

- [ ] Avoid unnecessary changes to prevent unnecessary updates across all pipelines.
- [ ] Verify if the template should skip code/files based on customisations (branding, modules).

## Command modifications

- [ ] Update all references to modified commands throughout documentation.
- [ ] Ensure mentions in nf-core/tools and nf-core/website documentation are current.
- [ ] Update the [`rich-codex.yml`](https://github.com/nf-core/website/blob/main/.github/rich-codex.yml) file in nf-core/website repository.
