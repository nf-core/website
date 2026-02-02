---
title: Reviewing nf-core/tools
subtitle: Review nf-core/tools pull requests
shortTitle: "Reviewing nf-core/tools"
markdownPlugin: checklist
---

nf-core/tools reviews ensure that changes to the tooling meet project standards and work reliably across the nf-core ecosystem.
When you review an nf-core/tools pull request, you examine the proposed changes and provide constructive feedback before maintainers merge them into the repository.

The infrastructure team oversees the review process for nf-core/tools, but community input helps catch issues and ensures changes work well in real-world scenarios.
Your perspective as a user of the tools is valuable, even if you're not familiar with all the internal code.

:::tip
Use this checklist for a quick reference whilst reviewing nf-core/tools PRs.
:::

## General checklist items

Check the fundamentals of the PR:

- [ ] Tests cover the added code (either old or newly added tests)
- [ ] Documentation exists or the PR adds it, including checking for corresponding PR updates to nf-core/website if needed
- [ ] Code does not hardcode "nf-core"; verify it uses `organisation` variable instead
- [ ] Code does not hardcode `modules`/`subworkflows`; verify it uses `component_type` instead

## Python styling checklist

Ensure the code follows Python best practices:

- [ ] Prefers `pathlib` over `os` module when possible for file operations
- [ ] Code uses typing
- [ ] Comments are sufficient and clear for readability
- [ ] Docstrings are updated (they generate API documentation automatically)
- [ ] If the PR adds a new command, verify that the API documentation skeleton includes it

## Template modifications

Template changes propagate to all pipelines through `nf-core pipelines sync`, so they require careful consideration:

- [ ] Avoid unnecessary changes to prevent unnecessary updates across all pipelines
- [ ] Verify if the template should skip code/files based on customisations (branding, modules)

## Command modifications

When the PR modifies commands, check that documentation is updated:

- [ ] All references to modified commands are updated throughout documentation
- [ ] Mentions in nf-core/tools and nf-core/website documentation are current
- [ ] The [`rich-codex.yml`](https://github.com/nf-core/website/blob/main/.github/rich-codex.yml) file in nf-core/website repository is updated

  :::note
  Command changes often require updates in multiple places.
  Missing documentation updates can confuse users who see different information in different locations.
  :::
