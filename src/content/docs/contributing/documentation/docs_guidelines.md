---
title: Writing style
subtitle: Guidelines for contributing to nf-core documentation.
weight: 10
parentWeight: 40
---

# nf-core/website - Documentation contribution guidelines

Many thanks for taking an interest in improving the nf-core documentation. To maintain consistency and clarity across our documentation, please follow these contribution guidelines.

## General writing guidelines

1. Use British English:

Write in British English, not US English. For example, "colour" instead of "color" and "organise" instead of "organize."

1. Use active voice:

Active voice makes sentences clearer and more direct.

Active: "nf-core curates pipelines."
Passive: "Pipelines are curated by nf-core."

1. Avoid gerunds:

Avoid using gerunds (verbs ending in -ing) to make headings and sentences more direct.

Example: "Run pipelines to produce results." not "Running pipelines produces results."

1. Be direct and brief:

Be direct and brief in your instructions and descriptions. Avoid words and qualifiers like _may_, _might_, _should_, _could_, _just_, or _even_.

Examples:

"To make changes to the pipeline", not "If you would like to make changes to the pipeline"
"Install with Conda", not "You can also install with Conda"

1. Avoid Latin abbreviations (e.g., i.e., etc.):

Write out the full phrase instead of using Latin abbreviations.

"For example," not "e.g."
"That is," not "i.e."

1. Avoid please and thank you:

Maintain a professional and direct tone by avoiding the use of "please" and "thank you".

## Formatting guidelines

1. Code and commands:

Use inline code formatting for commands, file names, and code snippets.

Example: Use `nextflow run nf-core/<pipeline_name>` to execute the pipeline.

1. Headings and titles:
   Use sentence case for headings and titles. Only the first word and proper nouns should be capitalized.

"Tuning workflow resources"
"Use a specific Nextflow version"

1. Lists:
   Use bullet points for lists where the order of items is not important. Use numbered lists for sequential steps.
   Use [checklists](/docs/contributing/website/markdown#checklist) if you're expecting someone to run through
   a set of instructions.

1. Hyperlinks:
   Write descriptive link text instead of using URLs directly.

:::tip
See [Markdown on the nf-core website](/docs/contributing/website/markdown#checklist) for instructions
on how to use markdown formatting and special elements such as "admonitions" for nice formatting.
You're reading text in a tip admonition right now!
:::

## Getting help

For further information or help, get in touch on the nf-core `docs` channel on [Slack](https://nf-co.re/join/slack/).
