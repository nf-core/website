---
title: Training material guide
subtitle: A guide to writing instructions for training material.
---

# Contributing a training guide to nf-core

Thank you for helping us write a guide to using an aspect of nf-core.

## Links

- [Questions](#questions)
- [Objectives](#objectives)
- [Glossary](#glossary)
- [Before I begin](#before-i-begin)
- [How to](#how-to)
  - [Bytesize](#bytesize-talk)
  - [Walkthrough](#walkthrough)
- [Key points](#key-points)
- [References](#references)

## Questions

- How can I contribute a training guide to nf-core?
- Where do I put the guide?
- How do I structure content?

## Objectives

- Learn how to structure an nf-core guide.
- Know where to submit the pull request.

## Glossary

<details>
    <summary>Click here to expand the glossary of potentially unfamiliar terms</summary>
    - **Learning Objective:** A predictive statement describing the specific outcomes that a training session is intended to achieve.
    - **Key point:** An important take home message.
    - **Cognitive load:** This relates to the amount of information that working memory can hold at one time.
    - **Fork a repository:** This makes a copy of the repository to your personal GitHub workspace. 
    - **Make a branch:** A branch is another version of the files in your repository, including a history of how files were changed from the start of the branch to the current state of files.
</details>

## Before I begin

A checklist of what you need to begin.

- Ask on #training about the guide you'd like to contribute.
- Make an issue on https://github.com/nf-core/nf-co.re/issues describing what you want to do, and assign yourself.
- GitHub account.

## How to

### Bytesize talk

There is no bytesize talk at present, but if you'd like present one, please make a suggestion on the [Slack #bytesize channel](https://nfcore.slack.com/channels/bytesize).

### Walkthrough

#### Use this as a template.

Use this file `https://github.com/nf-core/nf-co.re/markdown/developers/training_guide.md` as a template to guide you.

1. If you've not done so, make a fork of the [nf-core website repository](https://github.com/nf-core/nf-co.re/).
2. Make a new branch from the `master` branch in your fork named after the guide you indend to submit e.g., `gitpod_training_guide`.
3. Make a copy of the `training_guide.md` with an appropriate name to use as your skeleton.

#### Add questions.

Include between 1-3 questions your guide is intended to answer.

#### Add objectives.

Objectives are intended to describe the outcomes someone reading your guide should be able to do once they've read through it.
The objectives also function as a checklist of things that should be part of the guide.

#### Add a glossary.

The glossary should define potentially unfamilar terms to a novice trying to use your guide for the first time.

#### Prepare your reader.

Make a checklist of what your reader needs to have done before they begin using this guide.

#### Add the content.

Under the `### Guide` section, add a step by step walkthrough of what you intend to demonstrate. Remember to include visual aids if it will help understanding.

Initially it is unlikely there will be a bytesize talk to link to, so leave a message that one can be contributed.
If a bytesize talk exists, then provide a link using this html under the `### Bytesize talk` section. Update the `src` URL to the correct youtube link.

```
<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>
```

#### Add key points.

Key points serve as a summary of what to do, and are important take home messages from the guide.

#### Add supporting links.

Provide links to where the reader can get additional support, or useful further reading.

#### Format your document.

Format your document using the `prettier` VS code extension. This can be performed either in VS code, or in Gitpod.

## Key points

- Fork the nf-core website repository, and use the `training_guide.md` as a template.
- Add questions, objectives, and glossary terms.
- Have a checklist of what is needed before one begins.
- Provide a walkthrough, and a link to the bytesize talk if one exists.
- Provide a summary as key points.
- Provide links to support and documentation.
- Format your markdown document using the `prettier` VS code extension.

## References

Here are some additional links for where to find support or suggestions for further reading.

- [The `#training` Slack channel](https://nfcore.slack.com/channels/training)
- [The Turing Way Guide to collaboration](https://the-turing-way.netlify.app/collaboration/collaboration.html): A guide to effective collaboration.
- [Carpentries Handbook](https://docs.carpentries.org/): Describes The Carpentries way of training, including important concepts trainers should know.
- [Elixir train the trainer](https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT): Describes learning principles and how they apply to training.
