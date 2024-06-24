---
title: PR Review Checklist for nf-core/tools
subtitle: Guidelines for reviewing PRs in the nf-core/tools repository
shortTitle: nf-core/tools PR
markdownPlugin: checklist
---

While the [infrastructure team](https://github.com/orgs/nf-core/teams/infrastucture/members) is responsible for overseeing the PR review process for nf-core/tools, these guidelines can assist community members in reviewing PRs and ensure that the review process is consistent and effective. The following is a collection of points to have into account during the review process.

:::note
The nf-core/tools documentation is located in the nf-core/website. If the changes made to the nf-core/tools
affect this documentation, the contributor will have to open a PR to update it on that repo.
You will have to review that both updates are in sync.
:::

## General

- [ ] Check that the added code is covered by tests (either by old or newly added tests).
- [ ] Check that documentation exists or is added.
  - [ ] Check if a documentation PR to the [nf-core/website](https://github.com/nf-core/website) exists if needed.
- [ ] Check that `nf-core` is not hardcoded, use `organisation` instead.
- [ ] Check that `modules`/`subworkflows` is not hardcoded, use `component_type` instead.

## Python styling

- [ ] Check that [`pathlib`](https://docs.python.org/3/library/pathlib.html) is used instead of [`os`](https://docs.python.org/3/library/os.html) when possible.
- [ ] Check that the code uses [typing](https://docs.python.org/3/library/typing.html).
- [ ] Check that code comments are sufficient enough and clear.
- [ ] Check that docstrings are updated. This is specially important for linting, as the API documentation is automatically generated from these strings.
- [ ] If a new command was added, check that it is part of the API documentation skeleton, by running the `docs/api/generate-api-docs.sh` script

## Template modifications

If the PR is modifying the files under `nf-core/pipeline-template`:

- [ ] Check that unnecessary changes are NOT ADDED, to avoid template updates to all pipelines.
- [ ] Check if the code/file has to be skipped with any of the template customisations (branding, modules…)?

:::tip
You can see [an example of code](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/README.md?plain=1#L1-L10) which is skipped if the pipeline is not an nf-core pipeline.
And [an example of files](https://github.com/nf-core/tools/blob/master/nf_core/create.py#L61C9-L74) skipped depending of customisation options.
:::

## Command modifications

If the PR modifies how an nf-core command works.

- [ ] Check if all references to this command are updated.
- [ ] Check that mentions of this command on the nf-core/tools or nf-core/website documentation are updated.
- [ ] Make sure that the command in [rich-codex.yml](https://github.com/nf-core/website/blob/main/.github/rich-codex.yml) from the nf-core/website is updated.
