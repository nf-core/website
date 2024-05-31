---
title: Writing style
subtitle: Guidelines for contributing to nf-core documentation.
weight: 10
parentWeight: 40
markdownPlugin: addNumbersToHeadings
---

# Documentation contribution guidelines

Many thanks for taking an interest in improving the nf-core documentation. To maintain consistency and clarity across our documentation, please follow these contribution guidelines.

## General writing guidelines

### Use British English

Write in British English, not US English.

- ✅ _"colour"_ instead of ❌ _"color"_
- ✅ _"organise"_ instead of ❌ _"organize"_

### Use active voice

Active voice makes sentences clearer and more direct.

- ✅ Active: _"nf-core curates pipelines."_
- ❌ Passive: _"Pipelines are curated by nf-core."_

### Avoid gerunds:

Avoid using gerunds (verbs ending in -ing) to make headings and sentences more direct.

- ✅ _"Run pipelines to produce results."_
- ❌ _"Running pipelines produces results."_

### Be direct and brief

Be direct and brief in your instructions and descriptions. Avoid words and qualifiers like _may_, _might_, _should_, _could_, _just_, or _even_.

- ✅ _"To make changes to the pipeline"_
- ❌ _"If you would like to make changes to the pipeline"_
- ✅ _"Install with Conda"_
- ❌ _"You can also install with Conda"_

### Avoid Latin abbreviations (e.g., i.e., etc.)

Write out the full phrase instead of using Latin abbreviations.

- ✅ _"For example,"_ not ❌ _"e.g."_
- ✅ _"That is,"_ not ❌ _"i.e."_

### Avoid please and thank you

Maintain a professional and direct tone by avoiding the use of ❌ _"please"_ and ❌ _"thank you"_.

## Formatting guidelines

:::tip
See [Markdown on the nf-core website](/docs/contributing/website/markdown#checklist) for instructions
on how to use markdown formatting and special elements such as "admonitions" for nice formatting.
You're reading an admonition right now!
:::

### Code and commands

Use inline code formatting for commands, file names, and code snippets.

- ✅ _"Use `nextflow run nf-core/<pipeline_name>` to execute the pipeline."_

### Headings and titles

Use sentence case for headings and titles. Only the first word and proper nouns should be capitalized.

- ✅ _"Tuning workflow resources"_
- ✅ _"Use a specific Nextflow version"_

### Lists

- Use bullet points for lists where the order of items is not important.
- Use numbered lists for sequential steps.
- Use [numbered headings](/docs/contributing/documentation/markdown#numbered-headings) for guidelines where it is useful to be able to easily refer to a specific section of the docs (like this page!)
- Use [checklists](/docs/contributing/website/markdown#checklist) if you're expecting someone to run through
  a set of instructions.

### Hyperlinks

Write descriptive link text instead of using URLs directly.

- ✅ _"See the [nf-core website](https://nf-co.re)"_
- ❌ _"See the nf-core website: <https://nf-co.re>"_

## Getting help

For further information or help, get in touch on the nf-core [`#docs` channel](https://nfcore.slack.com/archives/C06RR8F5L3E) on [Slack](https://nf-co.re/join/slack/).
