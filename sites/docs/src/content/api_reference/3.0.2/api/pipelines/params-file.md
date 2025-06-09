# nf_core.params_file

Create a YAML parameter file

### _`class{:python}`_`nf_core.pipelines.params_file.ParamsFileBuilder(pipeline=None, revision=None){:python}`

Bases: `object`

Class to hold config option to launch a pipeline.

- **Parameters:**
  - **pipeline** (_str_ _,_ _optional_) – Path to a local pipeline path or a remote pipeline.
  - **revision** (_str_ _,_ _optional_) – Revision of the pipeline to use.

#### `format_group(definition, show_hidden=False){:python}`

Format a group of parameters of the schema as commented YAML.

- **Parameters:**
  - **definition** (_dict_) – Definition of the group from the schema
  - **show_hidden** (_bool_) – Whether to include hidden parameters
- **Returns:**
  Formatted output for a group
- **Return type:**
  str

#### `format_param(name, properties, required_properties=(), show_hidden=False){:python}`

Format a single parameter of the schema as commented YAML

- **Parameters:**
  - **name** (_str_) – Name of the parameter
  - **properties** (_dict_) – Properties of the parameter
  - **required_properties** (_list_) – List of required properties
  - **show_hidden** (_bool_) – Whether to include hidden parameters
- **Returns:**
  Section of a params-file.yml for given parameter
  None: If the parameter is skipped because it is hidden and show_hidden is not set
- **Return type:**
  str

#### `generate_params_file(show_hidden=False){:python}`

Generate the contents of a parameter template file.

Assumes the pipeline has been fetched (if remote) and the schema loaded.

- **Parameters:**
  **show_hidden** (_bool_) – Whether to include hidden parameters
- **Returns:**
  Formatted output for the pipeline schema
- **Return type:**
  str

#### `get_pipeline(){:python}`

Prompt the user for a pipeline name and get the schema

#### `write_params_file(output_fn='nf-params.yaml', show_hidden=False, force=False){:python}`

Build a template file for the pipeline schema.

- **Parameters:**
  - **output_fn** (_str_ _,_ _optional_) – Filename to write the template to.
  - **show_hidden** (_bool_ _,_ _optional_) – Include parameters marked as hidden in the output
  - **force** (_bool_ _,_ _optional_) – Whether to overwrite existing output file.
- **Returns:**
  True if the template was written successfully, False otherwise
- **Return type:**
  bool

### `nf_core.pipelines.params_file._print_wrapped(text, fill_char='-', mode='both', width=80, indent=0, drop_whitespace=True){:python}`

Helper function to format text for the params-file template.

- **Parameters:**
  - **text** (_str_) – Text to print
  - **fill_char** (_str_ _,_ _optional_) – Character to use for creating dividers. Defaults to ‘-‘.
  - **mode** (_str_ _,_ _optional_) – Where to place dividers. Defaults to “both”.
  - **width** (_int_ _,_ _optional_) – Maximum line-width of the output text. Defaults to 80.
  - **indent** (_int_ _,_ _optional_) – Number of spaces to indent the text. Defaults to 0.
  - **drop_whitespace** (_bool_ _,_ _optional_) – Whether to drop whitespace from the start and end of lines.
