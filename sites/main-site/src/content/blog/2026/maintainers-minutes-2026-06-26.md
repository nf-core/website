---
title: "Maintainers Minutes: May - June 2026"
subtitle: "Keeping you informed of the latest maintainers discussions"
pubDate: 2026-06-30T10:00:00+01:00
headerImage: "/assets/images/blog/maintainers-minutes-2024-07-26/maintainers-wide.png"
headerImageAlt: "Cartoon yellow rubber duck with nf-core logo badge on its body with the nf-core logo."
embedHeaderImage: false
authors:
  - "maxulysse"
label:
  - "maintainers"
---

The 'Maintainers Minutes' aims to give further insight into the workings of the [nf-core maintainers team](/governance#maintainers) by providing brief summaries of the monthly team meetings.

## Overview

The last two maintainer meetings covered a wide range of topics from AI related topics to infrastructure updates, and of course, the usual discussions about the meaning of life.

Major topics covered in this period:

- AI in general
- Configs strict syntax
- Deprecations documentation
- Topic migration progress
- Megatests snapshots
- Local modules in pipelines
- r-nf-core-utils
- Summer break and miscellaneous updates

### AI in general

#### RFC on agents

The [agents.md and AI steering document](https://github.com/nf-core/proposals/issues/141) RFC on agents has been accepted, and Igor is leading the effort.
This is a significant step towards formalizing how we handle AI-generated contributions and maintaining code quality standards.

#### Dealing with PRs made by/with agents

Related to that, we had a long discussion about how to deal with PRs made by or with agents.
These were some of the several key points that were raised, and what we talked about:

- Snapshots hallucination: AI agents can sometimes generate incorrect or hallucinated content
- Human responsibility: The human launching agents should be responsible, not the reviewers
- Agent origin doesn't matter: We don't care if contributions come from an agent or a human, but the volume is concerning
- Closing unsatisfactory PRs: Maxime recommends just closing PRs if the reviewer is not happy
- Assisted by tags: It would be nice to have an "assisted by" comment or tag somewhere to identify AI-generated contributions

The general consensus was that an "assisted by" tag would be useful, and that in the end the human launching the agent should be responsible for the quality of the contribution.
The reviewers should not be expected to verify the correctness of AI-generated content, and they are free to close PRs if they are not satisfied with the quality.

#### Need for more specifications and documentation

We discussed the need for more documentation and specifications about tools.
This could go in `meta.yml` files so that AI and humans don't remove or misunderstand the underlying logic of modules.

#### Igor is our only hope against AI

Igor will lead the effort for creating a steering file for agents, helping guide AI contributions to maintain our quality standards.

#### Open for testing AGENTS.md

AGENTS.md is progressing, and we can probably test it after summer.
Igor is looking for volunteers, and Maxime is volunteering to test it.

### Configs strict syntax

The configs repository is almost fully Nextflow lint valid with strict syntax.
This is a major milestone for improving code quality and consistency across nf-core configurations.

### Deprecations documentation

We need to reach a consensus on deprecation documentation.
Related issues include:

- [nf-core/modules#5602](https://github.com/nf-core/modules/issues/5602)
- [nf-core/tools#4260](https://github.com/nf-core/tools/issues/4260)

The consensus is to create a failure assertion in the test.
And we are assessing the modules deletion process a year after the deprecation.
nf-core tools should handle this nicely, as the modules.json file might still be populated with the deprecated module.

### Topic migration progress

End of May, the topic migration has reached the end of letter H, showing steady progress.
By the end of June, we are now at the letter M, which shows the dedication of Louis into leading this herculean task.

### Megatests snapshots

Maxime is hacking his way into megatests snapshots, working on improving nf-core's testing infrastructure to better handle large-scale testing scenarios.

### Local modules in pipelines

We discussed what to do with local modules in release PRs.
The consensus was to:

- Document them in an issue
- Try to push as many as possible to the nf-core/modules repo
- Recommend using `modules patch` when they diverge too much from the nf-core/modules repo

### r-nf-core-utils

The r-nf-core-utils package has been approved, and released.
This package will make writing R code for nf-core modules easier by providing helper functions for parsing Nextflow inputs and producing log/version files.
A proof of concept is available at [nf-core/modules#11933](https://github.com/nf-core/modules/pull/11933).

### Containers in nf-core

#### tools dev has a new container subcommand

Tools dev has a new subcommand for containers (in modules).
It will facilitate automated container building and management.
Modules will get ARM support in a much easier way.

#### Wave with different build template

Seqera containers via wave can be done using a different build template (pixi).
There might be some issues with the default base image that are being looked into.

### Meta pipelines

Meta pipelines are being discussed in [nextflow-io/nextflow#7213](https://github.com/nextflow-io/nextflow/pull/7213).

### Summer break

Summer break is coming, and both bytesizes and helpdesk will be on break.
The maintainers meeting is also under evaluation.

### End

That's all folks! We'll be back after the summer break with more updates and discussions.

\- :heart: from your #maintainers team!
