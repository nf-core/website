---
title: "Propose your pipeline before you build it"
subtitle: "Building in private and submitting a finished pipeline is the hardest way into nf-core"
headerImage: "https://images.unsplash.com/photo-1503387762-592deb58ef4e"
headerImageAlt: "A person drawing an architectural plan on paper with a ruler"
embedHeaderImage: true
pubDate: 2026-11-13T10:00:00+01:00
authors:
  - "ewels"
label:
  - "community"
  - "pipelines"
---

nf-core is quite an unusual community.
We build one set of best-practice pipelines that the community owns together, it is not simply a registry for listing compliant pipelines.
The people who maintain its pipelines share no employer and no funding body.
Thousands of modules get written once and reused across hundreds of workflows.
Pipelines outlive the PhD projects that started them.
Our most recent paper in [Genome Biology](https://doi.org/10.1186/s13059-025-03673-9) goes into this properly.
You can see it happening in older pipelines: the [contributor graph for nf-core/rnaseq](https://ewels.github.io/contributor-graphs/rnaseq.html) shows development passing between people over eight years as jobs and interests change.

To make this work, we have [guidelines](/docs/specifications/pipelines/overview) for how pipelines get built.
The [Specificity](/docs/specifications/pipelines/requirements/workflow_specificity) guideline says that nf-core carries one pipeline per data and analysis type, rather than having several pipelines that overlap.
[Community ownership](/docs/specifications/pipelines/requirements/community_owned) states that a pipeline belongs to the community rather than to the lab that started it, so when someone moves on, other people can pick it up.

It's far easier to satisfy the nf-core guidelines if you start the conversation early.
We encourage people to submit a pipeline proposal before starting any code, and then to build in public, together with the community.

## What building in private costs you

It's not unusual to have people proposing already-completed pipelines to nf-core.
Most proposals arrive with some code already written: 62% link to existing work, with varying degrees of completeness.
You can do that, but [building with the community](/docs/specifications/pipelines/recommendations/build_with_community) is a written recommendation, and it's central to how nf-core works.

If you build in private, you take on two risks.
The first is duplication.
nf-core allows [one pipeline per data and analysis type](/docs/specifications/pipelines/requirements/workflow_specificity), so if a pipeline with a similar purpose already exists or is proposed before yours, [you will be asked to contribute to that one instead, even if yours is more complete](/docs/contributing/contribute-new-pipelines#propose-before-you-build).

The second is feedback.
Arriving with a finished pipeline means that you might face fundamental design comments when acting on them costs the most, compared to having to include them at an earlier development stage.
There are fairly straightforward requirements such as [using the nf-core template](/docs/specifications/pipelines/requirements/use_the_template), [passing the lint tests](/docs/specifications/pipelines/requirements/linting) and being released under an [MIT licence](/docs/specifications/pipelines/requirements/mit_license), but code review itself can be more of an art than a science at times and often uncovers unforeseen issues or approaches which may need significant work.

These recommendations come from experience.
Since the [#new-pipelines](https://nfcore.slack.com/archives/CE6SDEDAA) Slack channel opened in November 2018, nf-core has decided on more than 300 pipeline proposals.
About three in five were accepted.
Of the roughly 120 that were not, about 100 were asked to join an existing pipeline instead.
Four in five of those had already written code.

<img src="/images/blog/propose-before-you-build/nf-core-proposals-all-time-light.svg" class="hide-dark" alt="Outcome of nf-core pipeline proposals, Nov 2018 to Nov 2026" />
<img src="/images/blog/propose-before-you-build/nf-core-proposals-all-time-dark.svg" class="hide-light" alt="Outcome of nf-core pipeline proposals, Nov 2018 to Nov 2026" />

## What works instead

Propose the idea while it is still an idea.
Work out the scope with people who know the field, including whether it is [the right size](/docs/specifications/pipelines/requirements/workflow_size) to be a pipeline at all.
Contribute your modules to [nf-core/modules](https://github.com/nf-core/modules) one at a time, each reviewed as it goes in, then assemble them into a pipeline.
Avoid building a pipeline mostly from local modules that nobody else can reuse.
By the time you have something running end to end you have had a long run of small reviews instead of one enormous one, and you have usually picked up a few people who care whether it works.

## "Does the idea stop being mine?"

This approach can feel quite alien to many.
It's understandable to worry that you may lose control of the development and end up with a pipeline that doesn't fit the needs of the group.

Shared ownership does not necessitate lack of control.
As maintainers, you decide where the pipeline goes.
In practice you will likely be the only developer, or at least the main one, with other people dropping in to add a small feature or fix a bug.

The [community owned](/docs/specifications/pipelines/requirements/community_owned) requirement sets out what you are agreeing to.
Other people may add features and make releases, and the community settles decisions by consensus.
The maintainers guide those decisions, with [final say from the core team and steering committee](/governance) if it comes to that.
In eight years and more than 150 pipelines, this has only been needed once or twice.
You are encouraged to develop using [a personal fork](/docs/specifications/pipelines/requirements/identity_branding) for your own work, as long as the nf-core repository stays the primary source.

Community ownership means development is not exclusive to one person or group, even when one person does all the maintenance.
If you move on and stop maintaining the pipeline, other people can pick it up and keep it working.

## Publishing in a journal

We actively want people to publish manuscripts about nf-core pipelines.
You can propose a pipeline in public before the paper or thesis about it comes out.
The primary developers of a pipeline are very welcome to publish it, and plenty already have, as the [publications list](/publications) shows.

The one thing that we do recommend is to talk about authorship early, if the pipeline has multiple contributors.
Say in the pipeline's Slack channel that you are writing a paper, so contributors know and you settle co-authorship before submission instead of after.
The [publication credit](/docs/specifications/pipelines/recommendations/publication_credit) recommendation and the [ownership and publication](/docs/contributing/contribute-new-pipelines#ownership-and-publication) docs cover the details, including how to cite nf-core and how to [mint a DOI](/docs/specifications/pipelines/recommendations/dois) for each release.

## If development really has to stay closed

Sometimes there are unavoidable reasons for pipeline development to be closed.
Not every pipeline built with nf-core tooling has to become an nf-core pipeline: you can [use the template, tools and components outside nf-core](/docs/developing/external-use) and propose later, knowing the duplication risk you are taking.
A turned-down proposal says nothing about the quality of your code or your science.
It is about what the community can realistically maintain.

## What we changed in the docs

This conversation comes up a lot, and it's difficult to convey the message to newcomers to the community.
To try to make it clearer, we have recently reorganised the documentation a little: [Contributing a new pipeline](/docs/contributing/contribute-new-pipelines) now starts with the proposal instead of with creating a repository, and carries a new [Propose before you build](/docs/contributing/contribute-new-pipelines#propose-before-you-build) section covering why to propose early, what community ownership means for your work, and what happens when you bring in a pipeline built elsewhere.
The warning now appears much earlier too, on the [contributing overview](/docs/contributing/overview), [how to contribute](/docs/contributing/how-to-contribute), the [what is nf-core](/docs/get_started/nf-core) page, and the specification pages people tend to land on first.

If you are weighing up a new pipeline, read that section and then come and talk to us.
Post in [Slack](/join), open a proposal, and ask the awkward questions about scope and licensing while they are still cheap to answer.
If a supervisor needs convincing, the paper and the contributor graphs above are a reasonable place to start.
Much better to have that conversation now, instead of after six months of work.

Welcome to nf-core.
We look forward to reading your pipeline proposals!
