# schema_lint

#### PipelineLint.schema_lint()

Pipeline schema syntax

Pipelines should have a `nextflow_schema.json` file that describes the different
pipeline parameters (eg. `params.something`, `--something`).

The lint test checks the schema for the following:

- Schema should be a valid JSON file
- Schema should adhere to [JSONSchema](https://json-schema.org/), Draft 7.
- Parameters can be described in two places:
  > - As `properties` in the top-level schema object
  > - As `properties` within subschemas listed in a top-level `definitions` objects
- The schema must describe at least one parameter
- There must be no duplicate parameter IDs across the schema and definition subschema
- All subschema in `definitions` must be referenced in the top-level `allOf` key
- The top-level `allOf` key must not describe any non-existent definitions
- Default parameters in the schema must be valid
- Core top-level schema attributes should exist and be set as follows:
  > - `$schema`: `https://json-schema.org/draft-07/schema`
  > - `$id`: URL to the raw schema file, eg. `https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json`
  > - `title`: `YOURPIPELINE pipeline parameters`
  > - `description`: The pipeline config `manifest.description`
- That the `input` property is defined and has a mimetype. A list of common mimetypes can be found [here](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types).

For example, an _extremely_ minimal schema could look like this:

```json
{
  "$schema": "https://json-schema.org/draft-07/schema",
  "$id": "https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json",
  "title": "YOURPIPELINE pipeline parameters",
  "description": "This pipeline is for testing",
  "properties": {
    "first_param": { "type": "string" }
  },
  "definitions": {
    "my_first_group": {
      "properties": {
        "second_param": { "type": "string" }
      }
    }
  },
  "allOf": [{ "$ref": "#/definitions/my_first_group" }]
}
```
