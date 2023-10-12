<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/params_file.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.params_file`

Create a YAML parameter file

## **Global Variables**

- **INTRO**
- **USAGE**
- **H1_SEPERATOR**
- **H2_SEPERATOR**

---

<a href="../../../../../../tools/nf_core/params_file.py#L78"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ParamsFileBuilder`

Class to hold config option to launch a pipeline.

**Args:**
pipeline (str, optional): Path to a local pipeline path or a remote pipeline. revision (str, optional): Revision of the pipeline to use.

<a href="../../../../../../tools/nf_core/params_file.py#L88"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline=None, revision=None)
```

Initialise the ParamFileBuilder class

**Args:**

- <b>`pipeline`</b> (str, optional): Path to a local pipeline path or a remote pipeline.
- <b>`revision`</b> (str, optional): Revision of the pipeline to use.

---

<a href="../../../../../../tools/nf_core/params_file.py#L135"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `format_group`

```python
format_group(definition, show_hidden=False)
```

Format a group of parameters of the schema as commented YAML.

**Args:**

- <b>`definition`</b> (dict): Definition of the group from the schema
- <b>`show_hidden`</b> (bool): Whether to include hidden parameters

**Returns:**

- <b>`str`</b>: Formatted output for a group

---

<a href="../../../../../../tools/nf_core/params_file.py#L174"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `format_param`

```python
format_param(name, properties, required_properties=(), show_hidden=False)
```

Format a single parameter of the schema as commented YAML

**Args:**

- <b>`name`</b> (str): Name of the parameter
- <b>`properties`</b> (dict): Properties of the parameter
- <b>`required_properties`</b> (list): List of required properties
- <b>`show_hidden`</b> (bool): Whether to include hidden parameters

**Returns:**

- <b>`str`</b>: Section of a params-file.yml for given parameter None: If the parameter is skipped because it is hidden and show_hidden is not set

---

<a href="../../../../../../tools/nf_core/params_file.py#L215"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `generate_params_file`

```python
generate_params_file(show_hidden=False)
```

Generate the contents of a parameter template file.

Assumes the pipeline has been fetched (if remote) and the schema loaded.

**Args:**

- <b>`show_hidden`</b> (bool): Whether to include hidden parameters

**Returns:**

- <b>`str`</b>: Formatted output for the pipeline schema

---

<a href="../../../../../../tools/nf_core/params_file.py#L107"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_pipeline`

```python
get_pipeline()
```

Prompt the user for a pipeline name and get the schema

---

<a href="../../../../../../tools/nf_core/params_file.py#L246"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `write_params_file`

```python
write_params_file(output_fn='nf-params.yaml', show_hidden=False, force=False)
```

Build a template file for the pipeline schema.

**Args:**

- <b>`output_fn`</b> (str, optional): Filename to write the template to. show_hidden (bool, optional): Include parameters marked as hidden in the output
- <b>`force`</b> (bool, optional): Whether to overwrite existing output file.

**Returns:**

- <b>`bool`</b>: True if the template was written successfully, False otherwise

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
