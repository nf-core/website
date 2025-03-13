---
title: Create a parameter file
subtitle: Create a skeleton / template parameter file for a pipeline
shortTitle: create-params-file
weight: 30
---

Sometimes it is easier to manually edit a parameter file than to use the web interface or interactive commandline wizard
provided by `nf-core pipelines launch`, for example when running a pipeline with many options on a remote server without a graphical interface.

You can create a parameter file with all parameters of a pipeline with the `nf-core pipelines create-params-file` command.
This file can then be passed to `nextflow` with the `-params-file` flag.

This command takes one argument - either the name of a nf-core pipeline which will be pulled automatically,
or the path to a directory containing a Nextflow pipeline _(can be any pipeline, doesn't have to be nf-core)_.

The generated YAML file contains all parameters set to the pipeline default value along with their description in comments.
This template can then be used by uncommenting and modifying the value of parameters you want to pass to a pipline run.

Hidden options are not included by default, but can be included using the `-x`/`--show-hidden` flag.
