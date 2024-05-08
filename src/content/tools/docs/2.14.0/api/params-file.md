# nf\_core.params\_file

Create a YAML parameter file

### *`class{:python}`*`nf_core.params_file.ParamsFileBuilder(pipeline=None, revision=None){:python}`

Bases: `object`

Class to hold config option to launch a pipeline.

* **Parameters:**
  * **pipeline** (*str* *,* *optional*) – Path to a local pipeline path or a remote pipeline.
  * **revision** (*str* *,* *optional*) – Revision of the pipeline to use.

#### `format_group(definition, show_hidden=False){:python}`

Format a group of parameters of the schema as commented YAML.

* **Parameters:**
  * **definition** (*dict*) – Definition of the group from the schema
  * **show\_hidden** (*bool*) – Whether to include hidden parameters
* **Returns:**
  Formatted output for a group
* **Return type:**
  str

#### `format_param(name, properties, required_properties=(), show_hidden=False){:python}`

Format a single parameter of the schema as commented YAML

* **Parameters:**
  * **name** (*str*) – Name of the parameter
  * **properties** (*dict*) – Properties of the parameter
  * **required\_properties** (*list*) – List of required properties
  * **show\_hidden** (*bool*) – Whether to include hidden parameters
* **Returns:**
  Section of a params-file.yml for given parameter
  None: If the parameter is skipped because it is hidden and show\_hidden is not set
* **Return type:**
  str

#### `generate_params_file(show_hidden=False){:python}`

Generate the contents of a parameter template file.

Assumes the pipeline has been fetched (if remote) and the schema loaded.

* **Parameters:**
  **show\_hidden** (*bool*) – Whether to include hidden parameters
* **Returns:**
  Formatted output for the pipeline schema
* **Return type:**
  str

#### `get_pipeline(){:python}`

Prompt the user for a pipeline name and get the schema

#### `write_params_file(output_fn='nf-params.yaml', show_hidden=False, force=False){:python}`

Build a template file for the pipeline schema.

* **Parameters:**
  * **output\_fn** (*str* *,* *optional*) – Filename to write the template to.
  * **show\_hidden** (*bool* *,* *optional*) – Include parameters marked as hidden in the output
  * **force** (*bool* *,* *optional*) – Whether to overwrite existing output file.
* **Returns:**
  True if the template was written successfully, False otherwise
* **Return type:**
  bool

### `nf_core.params_file._print_wrapped(text, fill_char='-', mode='both', width=80, indent=0, drop_whitespace=True){:python}`

Helper function to format text for the params-file template.

* **Parameters:**
  * **text** (*str*) – Text to print
  * **fill\_char** (*str* *,* *optional*) – Character to use for creating dividers. Defaults to ‘-‘.
  * **mode** (*str* *,* *optional*) – Where to place dividers. Defaults to “both”.
  * **width** (*int* *,* *optional*) – Maximum line-width of the output text. Defaults to 80.
  * **indent** (*int* *,* *optional*) – Number of spaces to indent the text. Defaults to 0.
  * **drop\_whitespace** (*bool* *,* *optional*) – Whether to drop whitespace from the start and end of lines.
