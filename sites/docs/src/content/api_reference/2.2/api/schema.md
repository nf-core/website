# nf_core.schema

Code to deal with pipeline JSON Schema

### _`class{:python}`_`nf_core.schema.PipelineSchema{:python}`

Bases: `object`

Class to generate a schema object with
functions to handle pipeline JSON Schema

#### `add_schema_found_configs(){:python}`

Add anything that’s found in the Nextflow params that’s missing in the pipeline schema

#### `build_schema(pipeline_dir, no_prompts, web_only, url){:python}`

Interactively build a new pipeline schema for a pipeline

#### `build_schema_param(p_val){:python}`

Build a pipeline schema dictionary for an param interactively

#### `get_schema_defaults(){:python}`

Generate set of default input parameters from schema.

Saves defaults to self.schema_defaults
Returns count of how many parameters were found (with or without a default value)

#### `get_schema_path(path, local_only=False, revision=None){:python}`

Given a pipeline name, directory, or path, set self.schema_filename

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

#### `prompt_remove_schema_notfound_config(p_key){:python}`

Check if a given key is found in the nextflow config params and prompt to remove it if note

Returns True if it should be removed, False if not.

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

#### `save_schema(){:python}`

Save a pipeline schema to a file

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
Checks that the schema “$id”, “title” and “description” attributes match the piipeline config.
