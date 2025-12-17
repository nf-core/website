---
title: Specifications
subtitle: Specifications for writing nf-core components and pipelines
weight: 1
---

## Overview

The nf-core specifications define standards and best practices for developing robust, reproducible, and maintainable Nextflow components and pipelines. While these specifications are written primarily for nf-core contributors and are enforced in the nf-core modules and pipelines repositories, the principles and patterns documented here represent proven approaches to writing high-quality bioinformatics workflows that can benefit the broader Nextflow community.

These specifications are the result of collaborative development by the nf-core community. They address common challenges in workflow development, from managing software dependencies and computational resources to ensuring reproducibility and facilitating community contributions.

## Module Specifications

The following specifications define standards for developing nf-core modules:

- **[General](/developers/specifications/modules/general)** - Foundation for module development including input/output handling, `ext.args`, multi-tool piping, compression, version emission, and script templating.
- **[Naming conventions](/developers/specifications/modules/naming-conventions)** - Standards for naming module files, processes, parameters, functions, channels, and outputs.
- **[Input/output options](/developers/specifications/modules/input-output-options)** - Guidelines for defining input channels, output emissions, and handling optional inputs and outputs.
- **[Documentation](/developers/specifications/modules/documentation)** - Requirements for `meta.yaml` files including tool descriptions, keywords, and ontology integration.
- **[Module parameters](/developers/specifications/modules/module-parameters)** - Guidelines for parameter usage ensuring modules remain flexible and reusable across different pipeline contexts.
- **[Resource requirements](/developers/specifications/modules/resource-requirements)** - Standards for specifying computational resources through process labels and the `task` directive.
- **[Software requirements](/developers/specifications/modules/software-requirements)** - Guidelines for declaring software dependencies using Conda, Docker, and Singularity through BioContainers.
- **[Testing](/developers/specifications/modules/testing)** - Requirements for nf-test including snapshot testing, stub tests, and CI configuration.
- **[Miscellaneous](/developers/specifications/modules/misc)** - Code formatting standards including the "Harshil Alignment" format.
