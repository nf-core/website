<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/schema.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.schema`

Code to deal with pipeline JSON Schema

---

<a href="../../../../../../tools/nf_core/schema.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineSchema`

Class to generate a schema object with functions to handle pipeline JSON Schema

<a href="../../../../../../tools/nf_core/schema.py#L30"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

Initialise the object

---

<a href="../../../../../../tools/nf_core/schema.py#L453"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `add_schema_found_configs`

```python
add_schema_found_configs()
```

Add anything that's found in the Nextflow params that's missing in the pipeline schema

---

<a href="../../../../../../tools/nf_core/schema.py#L301"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_schema`

```python
build_schema(pipeline_dir, no_prompts, web_only, url)
```

Interactively build a new pipeline schema for a pipeline

---

<a href="../../../../../../tools/nf_core/schema.py#L480"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_schema_param`

```python
build_schema_param(p_val)
```

Build a pipeline schema dictionary for an param interactively

---

<a href="../../../../../../tools/nf_core/schema.py#L102"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_schema_defaults`

```python
get_schema_defaults()
```

Generate set of default input parameters from schema.

Saves defaults to self.schema_defaults Returns count of how many parameters were found (with or without a default value)

---

<a href="../../../../../../tools/nf_core/schema.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_schema_path`

```python
get_schema_path(path, local_only=False, revision=None)
```

Given a pipeline name, directory, or path, set self.schema_filename

---

<a href="../../../../../../tools/nf_core/schema.py#L541"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_web_builder_response`

```python
get_web_builder_response()
```

Given a URL for a Schema build response, recursively query it until results are ready. Once ready, validate Schema and write to disk.

---

<a href="../../../../../../tools/nf_core/schema.py#L365"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_wf_params`

```python
get_wf_params()
```

Load the pipeline parameter defaults using `nextflow config` Strip out only the params. values and ignore anything that is not a flat variable

---

<a href="../../../../../../tools/nf_core/schema.py#L510"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_web_builder`

```python
launch_web_builder()
```

Send pipeline schema to web builder and wait for response

---

<a href="../../../../../../tools/nf_core/schema.py#L131"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_input_params`

```python
load_input_params(params_path)
```

Load a given a path to a parameters file (JSON/YAML)

These should be input parameters used to run a pipeline with the Nextflow -params-file option.

---

<a href="../../../../../../tools/nf_core/schema.py#L77"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_lint_schema`

```python
load_lint_schema()
```

Load and lint a given schema to see if it looks valid

---

<a href="../../../../../../tools/nf_core/schema.py#L94"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_schema`

```python
load_schema()
```

Load a pipeline schema from a file

---

<a href="../../../../../../tools/nf_core/schema.py#L286"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_skeleton_schema`

```python
make_skeleton_schema()
```

Make a new pipeline schema from the template

---

<a href="../../../../../../tools/nf_core/schema.py#L436"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_remove_schema_notfound_config`

```python
prompt_remove_schema_notfound_config(p_key)
```

Check if a given key is found in the nextflow config params and prompt to remove it if note

Returns True if it should be removed, False if not.

---

<a href="../../../../../../tools/nf_core/schema.py#L396"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_schema_notfound_configs`

```python
remove_schema_notfound_configs()
```

Go through top-level schema and all definitions sub-schemas to remove anything that's not in the nextflow config.

---

<a href="../../../../../../tools/nf_core/schema.py#L410"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_schema_notfound_configs_single_schema`

```python
remove_schema_notfound_configs_single_schema(schema)
```

Go through a single schema / set of properties and strip out anything that's not in the nextflow config.

Takes: Schema or sub-schema with properties key Returns: Cleaned schema / sub-schema

---

<a href="../../../../../../tools/nf_core/schema.py#L122"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `save_schema`

```python
save_schema()
```

Save a pipeline schema to a file

---

<a href="../../../../../../tools/nf_core/schema.py#L172"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_default_params`

```python
validate_default_params()
```

Check that all default parameters in the schema are valid Ignores 'required' flag, as required parameters might have no defaults

---

<a href="../../../../../../tools/nf_core/schema.py#L158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_params`

```python
validate_params()
```

Check given parameters against a schema and validate

---

<a href="../../../../../../tools/nf_core/schema.py#L193"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_schema`

```python
validate_schema(schema=None)
```

Check that the Schema is valid

Returns: Number of parameters found

---

<a href="../../../../../../tools/nf_core/schema.py#L239"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_schema_title_description`

```python
validate_schema_title_description(schema=None)
```

Extra validation command for linting. Checks that the schema "$id", "title" and "description" attributes match the piipeline config.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
