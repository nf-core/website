# nf_core.pipelines.schema

Code to deal with pipeline JSON Schema

### _`class{:python}`_`nf_core.pipelines.schema.PipelineSchema{:python}`

Bases: `object`

Class to generate a schema object with
functions to handle pipeline JSON Schema

#### `_update_validation_plugin_from_config() → None{:python}`

#### `add_schema_found_configs(){:python}`

Add anything that’s found in the Nextflow params that’s missing in the pipeline schema
Update defaults if they have changed

#### `build_schema(pipeline_dir, no_prompts, web_only, url){:python}`

Interactively build a new pipeline schema for a pipeline

#### `build_schema_param(p_val){:python}`

Build a pipeline schema dictionary for an param interactively

#### `check_for_input_mimetype(){:python}`

Check that the input parameter has a mimetype

Common mime types: <https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types>

- **Returns:**
  The mimetype of the input parameter
- **Return type:**
  mimetype (str)
- **Raises:**
  **LookupError** – If the input parameter is not found or defined in the correct place

#### `del_schema_filename() → None{:python}`

#### `get_schema_defaults() → None{:python}`

Generate set of default input parameters from schema.

Saves defaults to self.schema_defaults
Returns count of how many parameters were found (with or without a default value)

#### `get_schema_filename() → str{:python}`

#### `get_schema_path(path: str | Path, local_only: bool = False, revision: str | None = None) → None{:python}`

Given a pipeline name, directory, or path, set self.schema_filename

#### `get_schema_types() → None{:python}`

Get a list of all parameter types in the schema

#### `get_web_builder_response(){:python}`

Given a URL for a Schema build response, recursively query it until results are ready.
Once ready, validate Schema and write to disk.

#### `get_wf_params(){:python}`

Load the pipeline parameter defaults using nextflow config
Strip out only the params. values and ignore anything that is not a flat variable

#### `launch_web_builder(){:python}`

Send pipeline schema to web builder and wait for response

#### `load_input_params(params_path){:python}`

Load a given a path to a parameters file (JSON/YAML)

These should be input parameters used to run a pipeline with
the Nextflow -params-file option.

#### `load_lint_schema(){:python}`

Load and lint a given schema to see if it looks valid

#### `load_schema(){:python}`

Load a pipeline schema from a file

#### `make_skeleton_schema(){:python}`

Make a new pipeline schema from the template

#### `markdown_param_table(properties, required, columns){:python}`

Creates a markdown table for params from jsonschema properties section

- **Parameters:**
  - **properties** (_dict_) – A jsonschema properties dictionary
  - **required** (_list_) – A list of the required fields.
    Should come from the same level of the jsonschema as properties
  - **columns** (_list_) – A list of columns to write
- **Returns:**
  A string with the markdown table
- **Return type:**
  str

#### `markdown_to_html(markdown_str){:python}`

Convert markdown to html

#### `print_documentation(output_fn=None, format='markdown', force=False, columns=None){:python}`

Prints documentation for the schema.

#### `prompt_remove_schema_notfound_config(p_key){:python}`

Check if a given key is found in the nextflow config params and prompt to remove it if note

Returns True if it should be removed, False if not.

#### `remove_schema_empty_definitions(){:python}`

Go through top-level schema remove definitions that don’t have
any property attributes

#### `remove_schema_notfound_configs(){:python}`

Go through top-level schema and all definitions sub-schemas to remove
anything that’s not in the nextflow config.

#### `remove_schema_notfound_configs_single_schema(schema){:python}`

Go through a single schema / set of properties and strip out
anything that’s not in the nextflow config.

Takes: Schema or sub-schema with properties key
Returns: Cleaned schema / sub-schema

#### `sanitise_param_default(param){:python}`

Given a param, ensure that the default value is the correct variable type

#### `save_schema(suppress_logging=False){:python}`

Save a pipeline schema to a file

#### _`property{:python}`_`schema_filename{:python}`_: str_

#### `schema_to_markdown(columns){:python}`

Creates documentation for the schema in Markdown format.

#### `set_schema_filename(schema: str) → None{:python}`

#### `validate_config_default_parameter(param, schema_param, config_default){:python}`

Assure that default parameters in the nextflow.config are correctly set
by comparing them to their type in the schema

#### `validate_default_params(){:python}`

Check that all default parameters in the schema are valid
Ignores ‘required’ flag, as required parameters might have no defaults

Additional check that all parameters have defaults in nextflow.config and that
these are valid and adhere to guidelines

#### `validate_params(){:python}`

Check given parameters against a schema and validate

#### `validate_schema(schema=None){:python}`

Check that the Schema is valid

Returns: Number of parameters found

#### `validate_schema_title_description(schema=None){:python}`

Extra validation command for linting.
Checks that the schema “$id”, “title” and “description” attributes match the pipeline config.

### `nf_core.pipelines.schema.strip_required(node){:python}`
