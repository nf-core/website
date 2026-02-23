---
title: Build your nf-core pipeline dashboard with Depictio
category: community
slack: https://nfcore.slack.com/archives/C0ACF0TPF5E
location: Heidelberg and online
image: ""
image_alt: ""
leaders:
  weber8thomas:
    name: Thomas Weber
    slack: https://nfcore.slack.com/team/weber8thomas
---

Build interactive analysis dashboards for your nf-core pipeline results — no web development required.

## Goal

Create ready-to-use dashboards that combine MultiQC reports, Plotly figures, metric cards, and interactive components for any nf-core pipeline. Depictio offers a no-code drag-and-drop interface, while also providing versionable YAML configuration files for reproducibility and sharing.

## What is Depictio?

[Depictio](https://depictio.github.io/depictio-docs/) is an open-source platform for building interactive data dashboards from pipeline results. It features drag-and-drop layout, light/dark modes, and export capabilities — designed for bioinformaticians who want to visualize and share results without writing frontend code.

- **Demo** (ampliseq dashboard): [demo.depictio.embl.org](https://demo.depictio.embl.org/)
- **Docs**: [depictio.github.io/depictio-docs](https://depictio.github.io/depictio-docs/)
- **Install locally** or use the hosted instance at [SciLifeLab Serve](https://serve.scilifelab.se/) (available to Swedish researchers and SciLifeLab collaborators)

## Why build a dashboard?

- **Aggregate MultiQC reports** across runs to track quality over time
- **Create shareable apps** for publications and collaborators
- **Replace one-off scripts** with reusable, interactive visualizations

## What participants will do

1. Choose an nf-core pipeline you work with (or use nf-core test/megatest datasets)
2. Write the YAML configuration that tells Depictio which pipeline outputs to register (file paths, formats, data types)
3. Design a dashboard combining relevant plots, tables, and metrics using the drag-and-drop editor
4. Share your dashboard template and configuration with the community

## Recommended preparation

- Familiarity with an nf-core pipeline and its output structure
- Optionally bring results from a pipeline run (or plan to use test/megatest datasets)
- Skim the [Depictio documentation](https://depictio.github.io/depictio-docs/)
- No frontend or web development experience needed

_We welcome contributors of all experience levels._
