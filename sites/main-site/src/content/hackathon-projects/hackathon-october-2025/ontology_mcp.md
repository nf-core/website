---
title: ontology MCP
category: special-interest-groups
intro_video: ""
slack: https://nfcore.slack.com/archives/C08FKRKPDFV
image: "/assets/images/events/2025/hackathon-barcelona/ontologies_hackathon_image.jpeg"
image_alt: "Cat wondering between two question marks"
leaders:
  mirpedrol:
    name: Júlia Mir Pedrol
    slack: "https://nfcore.slack.com/team/U02JS8EFMMJ"
---

## Project Aim

Since last year, nf-core modules have a field in the `meta.yml` file to describe file ontologies.
But deciding which ontology term corresponds to each file can be a tedious task.

This project aims to work on the development of an MCP with tools useful to assign ontology terms to tool and file description,
which will be useful to tag nf-core modules.

An MCP (Model Context Protocol) is a way to provide tools for AI agents to have contact with the world and get information.
The MCP developed in the scope of this project should be able to assign ontology terms based on the description of a tool or a file.
If a term is not available, it can suggest a new term to be added to EDAM.

Work has been started in collaboration with the EDAM ontology organisation and the Bioconductor community at the [EDAM organisation](https://github.com/edamontology/edammcp).

For more context, you can find a [mermaid diagram](https://github.com/edamontology/edammcp/issues/8) of the intended use of the MCP and how it relates with nf-core, and the [initial discussion](https://github.com/edamontology/edammcp/issues/2) planning the different developement steps.

## Goals

1. Work on implementing MCP tools for assigning EDAM ontology terms in module's `meta.yml`
