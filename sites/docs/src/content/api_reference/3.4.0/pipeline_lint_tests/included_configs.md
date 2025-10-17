# included_configs

#### `PipelineLint.included_configs(){:python}`

Check that the pipeline nextflow.config includes the pipeline custom configs.

If the include line is uncommented, the test passes.
If the include line is commented, the test fails.
If the include line is missing, the test warns.

Can be skipped by adding the following to the .nf-core.yml file:
lint:

> included_configs: False
