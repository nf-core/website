---
title: Components
subtitle: Specifications for writing nf-core components
weight: 1
parentWeight: 10
---

Components are the building blocks of nf-core pipelines. They consist of modules (wrappers around individual tools) and subworkflows (combinations of multiple modules that perform related tasks). These specifications define the standards for developing reusable, interoperable components that can be shared across the nf-core ecosystem.

While these specifications are mandatory for components contributed to the nf-core repository, they represent proven best practices for developing high-quality, maintainable Nextflow components that can benefit any workflow development project.

## Module specifications

The following specifications define standards for developing nf-core modules:

- **[General](/developers/specifications/modules/general):** Foundation for module development including input/output handling, `ext.args`, multi-tool piping, compression, version emission, and script templating.
- **[Naming conventions](/developers/specifications/modules/naming-conventions):** Standards for naming module files, processes, parameters, functions, channels, and outputs.
- **[Input/output options](/developers/specifications/modules/input-output-options):** Guidelines for defining input channels, output emissions, and handling optional inputs and outputs.
- **[Documentation](/developers/specifications/modules/documentation):** Requirements for `meta.yaml` files including tool descriptions, keywords, and ontology integration.
- **[Module parameters](/developers/specifications/modules/module-parameters):** Guidelines for parameter usage ensuring modules remain flexible and reusable across different pipeline contexts.
- **[Resource requirements](/developers/specifications/modules/resource-requirements):** Standards for specifying computational resources through process labels and the `task` directive.
- **[Software requirements](/developers/specifications/modules/software-requirements):** Guidelines for declaring software dependencies using Conda, Docker, and Singularity through BioContainers.
- **[Testing](/developers/specifications/modules/testing):** Requirements for nf-test including snapshot testing, stub tests, and CI configuration.
- **[Miscellaneous](/developers/specifications/modules/misc):** Code formatting standards including the "Harshil Alignment" format.

## Subworkflow specifications

The following specifications define standards for developing nf-core subworkflows:

- **[General](/developers/specifications/subworkflows/general):** Foundation for subworkflow development including minimum subworkflow size and version reporting channels.
- **[Naming conventions](/developers/specifications/subworkflows/naming-conventions):** Standards for naming subworkflow files, parameters, functions, channels, and input/output structures.
- **[Input/output options](/developers/specifications/subworkflows/input-output-options):** Guidelines for defining required input and output channels, and handling optional inputs.
- **[Subworkflow parameters](/developers/specifications/subworkflows/subworkflow-parameters):** Guidelines for parameter usage ensuring subworkflows remain flexible and reusable across different pipeline contexts.
- **[Documentation](/developers/specifications/subworkflows/documentation):** Requirements for documenting channel structures in code comments and `meta.yml` files.
- **[Testing](/developers/specifications/subworkflows/testing):** Requirements for nf-test including scope of testing, tags for dependent modules, assertions, and CI configuration.
- **[Miscellaneous](/developers/specifications/subworkflows/misc):** Code formatting standards including the "Harshil Alignment" format.
