---
title: Create a parameter file
subtitle: Create a skeleton / template parameter file for a pipeline
shortTitle: create-params-file
weight: 30
---

Manually editing a parameter file can be easier than using the web interface or interactive command line wizard provided by `nf-core pipelines launch`, for example when running a pipeline with many options on a remote server without a graphical interface.

Create a parameter file with all pipeline parameters using the `nf-core pipelines create-params-file` command.
Pass this file to `nextflow` with the `-params-file` flag.

The command takes one argument: either the name of an nf-core pipeline (which will be pulled automatically) or the path to a directory containing a Nextflow pipeline (can be any pipeline, not just nf-core).

The generated YAML file contains all parameters set to the pipeline default value with their description in comments.
Use this template by uncommenting and modifying the values of parameters you want to pass to a pipeline run.

Hidden options are not included by default.
Include them using the `-x` or `--show-hidden` flag.
