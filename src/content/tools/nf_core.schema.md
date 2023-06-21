<!-- markdownlint-disable -->

<a href="../../nf_core/schema.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.schema`
Code to deal with pipeline JSON Schema  



---

<a href="../../nf_core/schema.py#L27"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineSchema`
Class to generate a schema object with functions to handle pipeline JSON Schema 

<a href="../../nf_core/schema.py#L31"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

Initialise the object 




---

<a href="../../nf_core/schema.py#L762"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `add_schema_found_configs`

```python
add_schema_found_configs()
```

Add anything that's found in the Nextflow params that's missing in the pipeline schema 

---

<a href="../../nf_core/schema.py#L577"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_schema`

```python
build_schema(pipeline_dir, no_prompts, web_only, url)
```

Interactively build a new pipeline schema for a pipeline 

---

<a href="../../nf_core/schema.py#L788"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_schema_param`

```python
build_schema_param(p_val)
```

Build a pipeline schema dictionary for an param interactively 

---

<a href="../../nf_core/schema.py#L431"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_for_input_mimetype`

```python
check_for_input_mimetype()
```

Check that the input parameter has a mimetype 

Common mime types: https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types 



**Returns:**
 
 - <b>`mimetype`</b> (str):  The mimetype of the input parameter 



**Raises:**
 
 - <b>`LookupError`</b>:  If the input parameter is not found or defined in the correct place 

---

<a href="../../nf_core/schema.py#L148"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_schema_defaults`

```python
get_schema_defaults()
```

Generate set of default input parameters from schema. 

Saves defaults to self.schema_defaults Returns count of how many parameters were found (with or without a default value) 

---

<a href="../../nf_core/schema.py#L50"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_schema_path`

```python
get_schema_path(path, local_only=False, revision=None)
```

Given a pipeline name, directory, or path, set self.schema_filename 

---

<a href="../../nf_core/schema.py#L853"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_web_builder_response`

```python
get_web_builder_response()
```

Given a URL for a Schema build response, recursively query it until results are ready. Once ready, validate Schema and write to disk. 

---

<a href="../../nf_core/schema.py#L646"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_wf_params`

```python
get_wf_params()
```

Load the pipeline parameter defaults using `nextflow config` Strip out only the params. values and ignore anything that is not a flat variable 

---

<a href="../../nf_core/schema.py#L817"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_web_builder`

```python
launch_web_builder()
```

Send pipeline schema to web builder and wait for response 

---

<a href="../../nf_core/schema.py#L179"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_input_params`

```python
load_input_params(params_path)
```

Load a given a path to a parameters file (JSON/YAML) 

These should be input parameters used to run a pipeline with the Nextflow -params-file option. 

---

<a href="../../nf_core/schema.py#L79"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_lint_schema`

```python
load_lint_schema()
```

Load and lint a given schema to see if it looks valid 

---

<a href="../../nf_core/schema.py#L108"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_schema`

```python
load_schema()
```

Load a pipeline schema from a file 

---

<a href="../../nf_core/schema.py#L562"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_skeleton_schema`

```python
make_skeleton_schema()
```

Make a new pipeline schema from the template 

---

<a href="../../nf_core/schema.py#L520"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `markdown_param_table`

```python
markdown_param_table(properties, required, columns)
```

Creates a markdown table for params from jsonschema properties section 



**Args:**
 
 - <b>`properties`</b> (dict):  A jsonschema properties dictionary 
 - <b>`required`</b> (list):  A list of the required fields.  Should come from the same level of the jsonschema as properties 
 - <b>`columns`</b> (list):  A list of columns to write 



**Returns:**
 
 - <b>`str`</b>:  A string with the markdown table 

---

<a href="../../nf_core/schema.py#L556"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `markdown_to_html`

```python
markdown_to_html(markdown_str)
```

Convert markdown to html 

---

<a href="../../nf_core/schema.py#L458"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_documentation`

```python
print_documentation(
    output_fn=None,
    format='markdown',
    force=False,
    columns=None
)
```

Prints documentation for the schema. 

---

<a href="../../nf_core/schema.py#L745"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_remove_schema_notfound_config`

```python
prompt_remove_schema_notfound_config(p_key)
```

Check if a given key is found in the nextflow config params and prompt to remove it if note 

Returns True if it should be removed, False if not. 

---

<a href="../../nf_core/schema.py#L677"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_schema_empty_definitions`

```python
remove_schema_empty_definitions()
```

Go through top-level schema remove definitions that don't have any property attributes 

---

<a href="../../nf_core/schema.py#L704"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_schema_notfound_configs`

```python
remove_schema_notfound_configs()
```

Go through top-level schema and all definitions sub-schemas to remove anything that's not in the nextflow config. 

---

<a href="../../nf_core/schema.py#L719"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_schema_notfound_configs_single_schema`

```python
remove_schema_notfound_configs_single_schema(schema)
```

Go through a single schema / set of properties and strip out anything that's not in the nextflow config. 

Takes: Schema or sub-schema with properties key Returns: Cleaned schema / sub-schema 

---

<a href="../../nf_core/schema.py#L116"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `sanitise_param_default`

```python
sanitise_param_default(param)
```

Given a param, ensure that the default value is the correct variable type 

---

<a href="../../nf_core/schema.py#L170"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `save_schema`

```python
save_schema(suppress_logging=False)
```

Save a pipeline schema to a file 

---

<a href="../../nf_core/schema.py#L495"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `schema_to_markdown`

```python
schema_to_markdown(columns)
```

Creates documentation for the schema in Markdown format. 

---

<a href="../../nf_core/schema.py#L286"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_config_default_parameter`

```python
validate_config_default_parameter(param, schema_param, config_default)
```

Assure that default parameters in the nextflow.config are correctly set by comparing them to their type in the schema 

---

<a href="../../nf_core/schema.py#L220"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_default_params`

```python
validate_default_params()
```

Check that all default parameters in the schema are valid Ignores 'required' flag, as required parameters might have no defaults 

Additional check that all parameters have defaults in nextflow.config and that these are valid and adhere to guidelines 

---

<a href="../../nf_core/schema.py#L207"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_params`

```python
validate_params()
```

Check given parameters against a schema and validate 

---

<a href="../../nf_core/schema.py#L334"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_schema`

```python
validate_schema(schema=None)
```

Check that the Schema is valid 

Returns: Number of parameters found 

---

<a href="../../nf_core/schema.py#L382"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_schema_title_description`

```python
validate_schema_title_description(schema=None)
```

Extra validation command for linting. Checks that the schema "$id", "title" and "description" attributes match the piipeline config. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
