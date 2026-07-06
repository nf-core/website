---
title: "From Confused to Contributing: My First nf-core PR Journey"
subtitle: "What I learned adding DIAMOND-blastp to proteinannotator and what I wish I'd known beforehand"
pubDate: 2026-07-10T09:00:00.000+01:00
headerImage: "https://plus.unsplash.com/premium_vector-1730152496499-809f0f287d43?q=80&w=880&auto=format&fit=crop&ixlib=rb-4.1.0&ixid=M3wxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8fA%3D%3D"
headerImageAlt: "A person with a bunch of tangled hair for a head"
authors: ["tracelail"]
label: ["community"]
---

## Introduction

Before I submitted my first nf-core pull request, I thought nf-core was a repository for pipelines. You go to the website, find the workflow you need, download it, plug in your parameters, and run it. Like Ensembl but for analysis. Maybe there's a tutorial or two. That's it.

The moment things started to shift was when my mentor pointed me toward the nf-core Slack help desk. I was stuck on something early in my PR, and she told me to go ask there. Once I was in and started looking around, I realized there was an entire world behind the website: active channels, people debugging pipelines in real time, module requests, tooling discussions, hackathons, summits, and a genuinely active community. It wasn't a repository. It was an ecosystem.

Here's the story of what happened when I tried to add [DIAMOND](https://github.com/bbuchfink/diamond)-blastp to the [proteinannotator](https://nf-co.re/proteinannotator/dev/) pipeline without fully understanding what I was getting into: the nesting doll structure that surprised me, the separate tools I didn't know I needed, the mental models that finally made things click, and what I'd do differently if I started over today.

If you're a bioinformatics student thinking about making your first nf-core contribution, this is the post I wish I'd had.


## Part 1: What I Thought nf-core Was vs. What It Actually Is

**What I thought:** A website where you go to download pre-built pipelines. You find the RNA-seq pipeline, input your parameters, and run a clean analysis without piecing together DESeq2 yourself. Convenient. Curated. Done.

**What it actually is:** A layered ecosystem, and the path to understanding it is more like a spiral than a straight line.

Here's the thing about the Nextflow training courses: they're well designed. The struggle that most beginners have, including myself, is that they want to apply what they're learning to a real scenario right away. The material isn't built this way. The training does what a good teacher should do, teach you from the ground up. You learn to walk before you run. It's quality material but you need to be patient with it.

I was a few modules into the [Hello Nextflow training](https://training.nextflow.io/latest/hello_nextflow/) when I started my PR, jumping in before finishing, which in retrospect was my first mistake. The courses had already shown me how to write a process, how inputs and outputs connect, and how a basic module works. So when I needed [DIAMOND](https://github.com/bbuchfink/diamond) blastp in my pipeline, I did exactly what the training had prepared me to do: I started writing it from scratch.

What I didn't know was that I didn't have to start from scratch. There are more than 1,400 pre-built, community-maintained modules available to install with a single command. Someone had already written the DIAMOND blastp module. I just had to know that it existed, and how to get it. That knowledge was there in the [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/) advanced track, waiting for me once I got there. So I was writing code the hard way, using the foundational skills I'd just learned, completely unaware that the whole point of nf-core tools was to make that unnecessary.

That's not a flaw in the Nextflow training material. I was impatient and tried to run a real-world PR on beginner knowledge. I needed to fully understand the structure above the module level.

It works like a nesting doll:

* **Modules** sit at the base: single-tool wrappers, one process per tool
* **Subworkflows** group modules together into a logical unit
* **Pipelines** are built from subworkflows and modules into a complete workflow

The other thing I missed early on, entirely my own fault, was the tooling. nf-core tools (`pip install nf-core`) is its own CLI suite. nf-test is another separate install. Neither is part of Nextflow itself. The nf-core docs actually spell this out clearly with [step-by-step instructions](https://nf-co.re/docs/get_started/environment_setup/overview) for every tool, I just hadn't read that page yet.

The analogy that finally made the whole structure click for me came from Python:

In Python, when you write `import pandas as pd`, you're declaring a dependency. You don't write pandas, you just use it. Everything happens under the hood. In Nextflow, you're writing the pandas. You're defining how the tool runs, what it takes as input, what it produces, how it connects to the next step. And like Python imports, those definitions live at the top level, declared before anything executes.

``` shell
# The nf-core module is already written by the community, you just install it:
$ nf-core modules install diamond/blastp
```

``` groovy
# Just import the nf-core module to use it
include { DIAMOND_BLASTP } from '../../../modules/nf-core/diamond/blastp/main'

# The local module: For this one, no community version exists, so you write it yourself
# modules/local/diamondpreparetaxa.nf
process DIAMONDPREPARETAXA {
    input:
    val taxondmp_zip   // e.g. ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    output:
    path "taxa/nodes.dmp" , emit: taxonnodes
    path "taxa/names.dmp" , emit: taxonnames

    script:
    """
    mkdir -p taxa/
    wget -q ${taxondmp_zip}
    tar -xzf taxdump.tar.gz -C taxa/
    """
}
```

The difference sounds subtle but discovering it mid-PR is not.


## Part 2: The Journey

[My PR](https://github.com/nf-core/proteinannotator/pull/50) started where someone else had already done part of the work, and I didn't fully understand what I was looking at.

The proteinannotator pipeline had been part of a previous nf-core hackathon. I knew I was picking up on an existing PR, that part was clear. What I didn't realize was the specific state it was in, or what I was even looking at. There were two DIAMOND-related files already in the repository: `diamond/blastp` and `diamond/makedb`. At the time, I assumed someone had been writing them as part of the PR. I didn't yet know that these were pre-built nf-core community modules that had simply been installed with a single command.

What was actually needed: local modules that didn't exist yet, a subworkflow to orchestrate everything, integration into the existing pipeline, and tests for all of it.

At the same time I was doing the beginner [Hello Nextflow training](https://training.nextflow.io/latest/hello_nextflow/), working through them alongside the PR rather than before it. My first instinct when I opened the repository was to figure out where things went. I knew there were modules. I knew they nested somehow but the actual shape of it, `modules/nf-core/` versus `modules/local/`, the `subworkflows/` directory, the `nextflow.config`, the `nextflow_schema.json`, the `conf/` folder, took time to internalize.

There is a learning curve, but it's not insurmountable. Once the structure clicks, it clicks completely. You'll be stuck on something for days, and then suddenly it makes sense, and everything else falls into place around it. The key is not giving up during the steep part.

I am a test-driven developer by instinct. So my approach was: write a local module, write a test for it, get the test passing, then move on. That approach served me well. It meant I always had something concrete to check my progress against, and it forced me to understand what each process was actually producing before I tried to wire it into anything larger.

What it didn't protect me from was the layer confusion.

Nextflow has three distinct layers that can look similar on the surface. The `script:` block is where you write bash, Python, or whatever shell commands your tool needs, familiar stuff. Above that is the Nextflow language level, where channel operators like `.filter {}` and `.map {}` act on the dataflow/channels emitting values asynchronously to the next step in the pipeline. And woven through that are collection methods like `.findAll {}` and `.find {}`, which look almost identical to channel operators but act on in-memory data structures synchronously and return a value immediately rather than emitting to a channel.

For example, `.filter { item -> condition }` and `.findAll { item -> condition }` look nearly the same. Same closure syntax, similar semantics. But one acts on a data stream and one acts on a List, and mixing them up produces errors.

I was making changes on all three layers simultaneously, and when something went wrong, it was often difficult to tell which layer was causing the error I was looking at.

The advanced Nextflow training courses are what eventually untangled this. I assumed they were for people who had already submitted a few PRs or were writing major pipelines. They're not. They cover exactly this kind of architectural distinction: how channel operations work at the workflow level, how the script block is a separate runtime, why the schema and config need to match exactly. If I had done them earlier, this confusion would have been resolved sooner.

I had been pointed toward Slack by my mentor early on, specifically for the weekly help desk channel. The first time I really needed it, I was stuck on a test that kept failing with an error I couldn’t resolve. The issue turned out to be the meta map. In nf-core, almost every process input expects a tuple: a metadata map (typically containing an `id` field) plus the actual file or value. My test didn't include that meta map. I had just passed the file directly. Someone in the help desk was able to identify it in just a few minutes.

This is worth saying plainly: Stack Overflow and general searches won't reliably get you through Nextflow problems. It's a niche ecosystem and the current knowledge lives within the Seqera community. For me, the nf-core Slack weekly help desk channel was what got me through. It is staffed by volunteers who actually understand the framework and want to help. Every wall I hit that might have taken me hours to resolve was solved there in minutes. Other great places to look are the Seqera community forum at [community.seqera.io](https://community.seqera.io), where questions are searchable, and Seqera AI, which is trained on Nextflow and nf-core documentation.

The DIAMOND PR stayed untouched for nearly a year. When I came back for the final push recently, I had done my own self-guided Nextflow training in the interim. Returning with that foundation, the things that had previously confused me finally made sense. I finished the local modules, built the DIAMOND subworkflow that orchestrated everything, and integrated it into the existing `functional_annotation` subworkflow, which feeds into `main.nf`. I updated the schema and the config, wrote the documentation, ran the full lint suite, and submitted.

The review feedback was minimal in volume but highly directive in quality. Two rounds of review, concrete changes each time. The PR is now in its third and final review. It's not merged yet. But it's right there.


## Part 3: If I Started Over

Here's what I'd do differently:

* Browse a complete pipeline before you touch anything.
    * Do the Hello Nextflow training first, the Hello nf-core training afterwards, and then look at the nf-core/rnaseq pipeline on GitHub.
    * Spend fifteen minutes clicking through the directories. Don't try to understand it all. Just see the shape of a finished pipeline before you start building one.
* Set up your environment from day one.
    * Every nf-core module needs a container image to run.
    * Make sure Docker, Singularity, or Conda is working before you write a single line of code.
* Set up a dedicated environment for nf-core tools and nf-test immediately.
    * Both are separate installs from Nextflow itself.
* Do the advanced track Nextflow training courses early.
    * They are not meant for expert users only. They are for anyone who wants to understand how the system actually fits together.
* Go to Slack.
    * The weekly help desk channel is especially helpful for newcomers.
* Start with something manageable.
    * Look for a PR that's partially done, or scope your first contribution to adding a tool that already exists as an nf-core module. The configs, the testing, and the structure will be the real learning curve regardless of what you pick.


## Conclusion

When I started this PR, I thought nf-core was a repository where you went to download pipelines. What I learned is that it's an ecosystem and a community, and a genuinely welcoming one. Every time I got stuck, someone in the Slack community walked me through it. That's why I'm now a Nextflow ambassador, and why I'm writing this post.

You don't need to submit a PR to get started. You don't need to understand everything before you try. Run an existing pipeline. Go through training. See if it clicks. The community is there when you need it, and the framework is better than it looks from the outside.
