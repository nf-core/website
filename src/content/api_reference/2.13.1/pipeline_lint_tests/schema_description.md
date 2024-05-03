# schema_description

#### `PipelineLint.schema_description(){:python}`

Check that every parameter in the schema has a description.

The `nextflow_schema.json` pipeline schema should describe every flat parameter.

Furthermore warns about parameters outside of groups.

- Warning: Parameters in `nextflow_schema.json` without a description
- Warning: Parameters in `nextflow_schema.json` that are defined outside of a group
