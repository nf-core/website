# schema_params

#### `PipelineLint.schema_params(){:python}`

Check that the schema describes all flat params in the pipeline.

The `nextflow_schema.json` pipeline schema should describe every flat parameter
returned from the `nextflow config` command (params that are objects or more complex structures are ignored).

- Failure: If parameters are found in `nextflow_schema.json` that are not in `nextflow_schema.json`
- Warning: If parameters are found in `nextflow_schema.json` that are not in `nextflow_schema.json`
