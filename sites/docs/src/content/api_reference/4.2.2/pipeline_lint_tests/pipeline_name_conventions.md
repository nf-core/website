# pipeline_name_conventions

#### `PipelineLint.pipeline_name_conventions(){:python}`

Checks that the pipeline name adheres to nf-core conventions.

In order to ensure consistent naming, pipeline names should contain only lower case, alphanumeric characters.
Otherwise a warning is displayed.

:::warning
DockerHub is very picky about image names and doesn’t even allow hyphens (we are `nfcore`).
This is a large part of why we set this rule.
:::
