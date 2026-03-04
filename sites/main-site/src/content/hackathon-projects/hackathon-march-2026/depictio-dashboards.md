---
title: Build your nf-core pipeline dashboard with Depictio
category: community
slack: https://nfcore.slack.com/archives/C0AFS06G9HR
location: Online
image: "/assets/images/events/2026/hackathon-march/depictio_logo.png"
image_alt: "Depictio logo"
leaders:
  weber8thomas:
    name: Thomas Weber
    slack: https://nfcore.slack.com/team/weber8thomas
---

Build interactive analysis dashboards for your nf-core pipeline results — no web development required.

## Goal

Create ready-to-use dashboard templates that combine MultiQC reports, Plotly figures, metric cards, geomaps, images, and interactive cross-filtering for any nf-core pipeline. Once a template exists for a pipeline, anyone can populate it with their own results — no configuration needed.

## What is Depictio?

[Depictio](https://depictio.github.io/depictio-docs/) is an open-source platform for building interactive data dashboards from pipeline results. It features a drag-and-drop UI, YAML-based configuration, sharing capabilities, designed for bioinformaticians who want to visualize and share results without writing frontend code. 

:::note{title="See it in action"}
Check out the [ampliseq dashboard demo](https://demo.depictio.embl.org/).
:::

- **Docs**: [depictio.github.io/depictio-docs](https://depictio.github.io/depictio-docs/)
- **Install locally** or use the hosted instance at [SciLifeLab Serve](https://serve.scilifelab.se/) (available to Swedish researchers and SciLifeLab collaborators)

## Why build a dashboard?

- **Aggregate MultiQC reports** across runs to track quality over time
- **Create shareable apps** for publications and collaborators
- **Replace one-off scripts** with reusable, interactive visualizations

## Tasks

:::info{title="No coding required"}
This project is suitable for newcomers to nf-core and open source contributions. No frontend or web development experience needed.
:::

1. Choose an nf-core pipeline you love and/or work with
2. Write the YAML configuration that tells Depictio which pipeline outputs to register (file paths, formats, data types) and how they link together
3. Design a multi-tab dashboard combining relevant plots, tables, and metrics using the drag-and-drop editor
4. Share your dashboard template and configuration with the community

## Lightweight peer-review

To ensure template quality, we'd like domain experts to do a lightweight review: check whether the structure, selected metrics, and interpretation logic make sense from their perspective, without taking on any Depictio maintenance responsibility. A sanity check from someone who knows the biology, not a long-term commitment.

## Longer-term vision

The longer-term goal is to provide standardised, reusable templates that can be populated automatically via an nf-core plugin triggered at the end of a pipeline run, pushing workflow outputs directly into Depictio. Multiple templates may exist for the same pipeline depending on the type of analysis.

## Prerequisites

- Familiarity with an nf-core pipeline and its output structure (which files are produced, their formats) is the main prerequisite
- Optionally bring results from a pipeline run (or plan to use test/megatest datasets)
- Skim the [Depictio documentation](https://depictio.github.io/depictio-docs/)
- No frontend or web development experience needed. Familiarity with YAML is recommended.

_We welcome contributors of all experience levels._
