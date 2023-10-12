<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/schema_lint.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.schema_lint`

---

<a href="../../../../../../tools/nf_core/lint/schema_lint.py#L8"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `schema_lint`

```python
schema_lint()
```

Pipeline schema syntax

Pipelines should have a `nextflow_schema.json` file that describes the different pipeline parameters (eg. `params.something`, `--something`).

.. tip:: Reminder: you should generally never need to edit this JSON file by hand. The `nf-core schema build` command can create _and edit_ the file for you to keep it up to date, with a friendly user-interface for customisation.

The lint test checks the schema for the following:

- Schema should be a valid JSON file _ Schema should adhere to `JSONSchema <https://json-schema.org/>`\_, Draft 7. _ Parameters can be described in two places:

- As `properties` in the top-level schema object \* As `properties` within subschemas listed in a top-level `definitions` objects

- The schema must describe at least one parameter _ There must be no duplicate parameter IDs across the schema and definition subschema _ All subschema in `definitions` must be referenced in the top-level `allOf` key _ The top-level `allOf` key must not describe any non-existent definitions _ Default parameters in the schema must be valid \* Core top-level schema attributes should exist and be set as follows:

- `$schema`: `https://json-schema.org/draft-07/schema` _ `$id`: URL to the raw schema file, eg. `https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json` _ `title`: `YOURPIPELINE pipeline parameters` _ `description`: The pipeline config `manifest.description` _ That the `input` property is defined and has a mimetype. A list of common mimetypes can be found `here <https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types>`\_.

For example, an _extremely_ minimal schema could look like this:

.. code-block:: json

{ "$schema": "https://json-schema.org/draft-07/schema",  "$id": "https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json", "title": "YOURPIPELINE pipeline parameters", "description": "This pipeline is for testing", "properties": { "first_param": { "type": "string" } }, "definitions": { "my_first_group": { "properties": { "second_param": { "type": "string" } } } }, "allOf": [{"$ref": "#/definitions/my_first_group"}] }

.. tip:: You can check your pipeline schema without having to run the entire pipeline lint by running `nf-core schema lint` instead of `nf-core lint`

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
