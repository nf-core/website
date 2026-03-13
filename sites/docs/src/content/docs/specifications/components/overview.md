---
title: Components
subtitle: Write nf-core components
weight: 1
parentWeight: 10
---

Components are the building blocks of nf-core pipelines.
They consist of modules (wrappers around individual tools) and subworkflows (combinations of multiple modules that perform related tasks).
These specifications define the standards for developing reusable, interoperable components that can be shared across the nf-core ecosystem.

While these specifications are mandatory for components contributed to the nf-core repository, they represent proven best practices for developing high-quality, maintainable Nextflow components that can benefit any workflow development project.

## Module specifications

The following specifications define standards for developing nf-core modules:
- **[General](/docs/specifications/components/modules/general):** Foundation for module development including input/output handling, `ext.args`, multi-tool piping, compression, version emission, and script templating.
- **[Naming conventions](/docs/specifications/components/modules/naming-conventions):** Standards for naming module files, processes, parameters, functions, channels, and outputs.
- **[Input/output options](/docs/specifications/components/modules/input-output-options):** Guidelines for defining input channels, output emissions, and handling optional inputs and outputs.
- **[Documentation](/docs/specifications/components/modules/documentation):** Requirements for `meta.yaml` files including tool descriptions, keywords, and ontology integration.
- **[Module parameters](/docs/specifications/components/modules/module-parameters):** Guidelines for parameter usage ensuring modules remain flexible and reusable across different pipeline contexts.
- **[Resource requirements](/docs/specifications/components/modules/resource-requirements):** Standards for specifying computational resources through process labels and the `task` directive.
- **[Software requirements](/docs/specifications/components/modules/software-requirements):** Guidelines for declaring software dependencies using Conda, Docker, and Singularity through BioContainers.
- **[Testing](/docs/specifications/components/modules/testing):** Requirements for nf-test including snapshot testing, stub tests, and CI configuration.

## Subworkflow specifications

The following specifications define standards for developing nf-core subworkflows:

- **[General](/docs/specifications/components/subworkflows/general):** Foundation for subworkflow development including minimum subworkflow size and version reporting channels.
- **[Naming conventions](/docs/specifications/components/subworkflows/naming-conventions):** Standards for naming subworkflow files, parameters, functions, channels, and input/output structures.
- **[Input/output options](/docs/specifications/components/subworkflows/input-output-options):** Guidelines for defining required input and output channels, and handling optional inputs.
- **[Subworkflow parameters](/docs/specifications/components/subworkflows/subworkflow-parameters):** Guidelines for parameter usage ensuring subworkflows remain flexible and reusable across different pipeline contexts.
- **[Documentation](/docs/specifications/components/subworkflows/documentation):** Requirements for documenting channel structures in code comments and `meta.yml` files.
- **[Testing](/docs/specifications/components/subworkflows/testing):** Requirements for nf-test including scope of testing, tags for dependent modules, assertions, and CI configuration.
