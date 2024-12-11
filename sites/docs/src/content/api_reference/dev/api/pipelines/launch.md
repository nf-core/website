# nf_core.pipelines.launch

Launch a pipeline, interactively collecting params

### _`class{:python}`_`nf_core.pipelines.launch.Launch(pipeline=None, revision=None, command_only=False, params_in=None, params_out=None, save_all=False, show_hidden=False, url=None, web_id=None){:python}`

Bases: `object`

Class to hold config option to launch a pipeline

#### `build_command(){:python}`

Build the nextflow run command based on what we know

#### `get_pipeline_schema(){:python}`

Load and validate the schema from the supplied pipeline

#### `get_web_launch_response(){:python}`

Given a URL for a web-gui launch response, recursively query it until results are ready.

#### `launch_pipeline(){:python}`

#### `launch_web_gui(){:python}`

Send schema to nf-core website and launch input GUI

#### `launch_workflow(){:python}`

Launch nextflow if required

#### `merge_nxf_flag_schema(){:python}`

Take the Nextflow flag schema and merge it with the pipeline schema

#### `print_param_header(param_id, param_obj, is_group=False){:python}`

#### `prompt_group(group_id, group_obj){:python}`

Prompt for edits to a group of parameters (subschema in ‘definitions’)

- **Parameters:**
  - **group_id** – Parameter ID (string)
  - **group_obj** – JSON Schema keys (dict)
- **Returns:**
  val answers
- **Return type:**
  Dict of param_id

#### `prompt_param(param_id, param_obj, is_required, answers){:python}`

Prompt for a single parameter

#### `prompt_schema(){:python}`

Go through the pipeline schema and prompt user to change defaults

#### `prompt_web_gui(){:python}`

Ask whether to use the web-based or cli wizard to collect params

#### `sanitise_web_response(){:python}`

The web builder returns everything as strings.
Use the functions defined in the cli wizard to convert to the correct types.

#### `set_schema_inputs(){:python}`

Take the loaded schema and set the defaults as the input parameters
If a nf_params.json file is supplied, apply these over the top

#### `single_param_to_questionary(param_id, param_obj, answers=None, print_help=True){:python}`

Convert a JSONSchema param to a Questionary question

- **Parameters:**
  - **param_id** – Parameter ID (string)
  - **param_obj** – JSON Schema keys (dict)
  - **answers** – Optional preexisting answers (dict)
  - **print_help** – If description and help_text should be printed (bool)
- **Returns:**
  Single Questionary dict, to be appended to questions list

#### `strip_default_params(){:python}`

Strip parameters if they have not changed from the default
