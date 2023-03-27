---
title: Training material guide
subtitle: A guide to writing instructions for training material.
---

# Contributing a training guide to nf-core

Thank you for helping us write a guide to using an aspect of nf-core.

## Scope

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
  - **Summative assessment:** A way of evaluating the learners progress at the end of a lesson.
  - **Fork a repository:** This makes a copy of the repository to your personal GitHub workspace.
  - **Make a branch:** A branch is another version of the files in your repository, including a history of how files were changed from the start of the branch to the current state of files.

</details>

## Before I begin

A checklist of what you need to begin.

- Be familiar with Markdown. If you're new to markdown, have a go a this [Markdown tutorial](https://www.markdowntutorial.com/).
- Ask on the nf-core Slack [`#training` channel](https://nfcore.slack.com/channels/training) about the guide you'd like to contribute.
- Make an issue on [nf-core website issues](https://github.com/nf-core/nf-co.re/issues) describing what you want to do, and assign yourself.

## How to

### Bytesize talk

There is no bytesize talk about how to create training materials at present, but if you'd like present one, please make a suggestion on the Slack [`#bytesize` channel](https://nfcore.slack.com/channels/bytesize).

### Walkthrough

#### Use this as a template.

Use this file `https://github.com/nf-core/nf-co.re/markdown/docs/contributing/training_guide.md` as a template to guide you.

1. Make a fork of the [nf-core website repository](https://github.com/nf-core/nf-co.re/).
2. Make a new branch from the `master` branch in your fork and check it out. Name it after the guide you intend to submit. For example: `git checkout -b gitpod_training_guide`
3. Make a copy of the `markdown/docs/contributing/training_guide.md` with an appropriate name to use as your skeleton.

#### Add framing questions.

Include between 1-3 framing questions your guide is intended to answer under the `## Scope` header above.
Framing questions are open ended questions, that set the scope for the guide.

#### Add objectives.

Objectives are intended to describe the intended learning outcomes someone reading your guide should be able to do once they've read through it.
The objectives also function as a checklist of things that should be part of the guide.

#### Add a glossary.

The glossary should define potentially unfamiliar terms to a novice trying to use your guide for the first time.

#### Prepare your reader.

Make a checklist of what your reader needs to have done before they begin using this guide.

#### Add the content.

Under the `### Guide` section, add a step by step walkthrough of what you intend to demonstrate. Remember to include visual aids if it will help understanding.

Initially it is unlikely there will be a bytesize talk to link to, so leave a message that one can be contributed.
If a bytesize talk exists, then provide a link using this html under the `### Bytesize talk` section. Update the `src` URL to the correct youtube link.

```html
<div class="ratio ratio-16x9">
  <iframe
    width="560"
    height="315"
    src="https://www.youtube.com/embed/xuNYATGFuw4"
    title="YouTube video player"
    frameborder="0"
    allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
    allowfullscreen
  ></iframe>
</div>
```

#### Add key points.

Key points serve as a summary of what to do, and are important take home messages from the guide.

#### Add supporting links.

Provide links to where the reader can get additional support, or useful further reading.

#### Add how to check your progress

Guides should provide some form of summative assessment at the end if appropriate.
Provide 3-4 questions you believe would help the reader to check their understanding or learning progress, after they've studied the material.
Add either a multiple choice or a true/false type of question, and feedback in case of correct or incorrect answer.
You can use the HTML `<details>` and `<summary>` tags to hide answers, e.g.,

```html
Your question here.
<details>
  <summary>Click to expand the solution.</summary>

  Text of your solution here.
</details>
```

#### Format your document.

Markdown must be formatted according to the rules set by a linting tool called _Prettier_ (<https://prettier.io/>).
We recommend using an extension for your code editor to fix markdown for you automatically when you save a file - there are Prettier plugins [for most editors](https://prettier.io/docs/en/editors.html), including [VSCode](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode).
Alternatively you can install and run the Prettier command line tool: `prettier -w .`

#### Adding your guide to the website navigation

For now, new guides need to be added to the website navigation manually.
A member of the core team will do that for you.

#### Submit a pull request.

Follow the instructions on the Github Documentation for [making a pull request from a fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).
Link the issue created for the guide you wanted to submit, by adding it under the developers box on the right side (below the milestone box, and above the notifications box)

## Key points

- Fork the nf-core website repository, and use the `training_guide.md` as a template.
- Add scope, objectives, and glossary terms.
- Have a checklist of what is needed before one begins.
- Provide a walkthrough, and a link to the bytesize talk if one exists.
- Provide a summary as key points.
- Provide links to support and documentation.
- Provide questions to evaluate if the learner has understood the guide content if applicable.
- Format your markdown document using the `prettier` VS code extension.

## References

Here are some additional links for where to find support or suggestions for further reading.

- [The `#training` Slack channel](https://nfcore.slack.com/channels/training)
- [The Turing Way Guide to collaboration](https://the-turing-way.netlify.app/collaboration/collaboration.html): A guide to effective collaboration.
- [Carpentries Handbook](https://docs.carpentries.org/): Describes The Carpentries way of training, including important concepts trainers should know.
- [Elixir train the trainer](https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT): Describes learning principles and how they apply to training.
