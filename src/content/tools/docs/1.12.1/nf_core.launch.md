<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/launch.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.launch`

Launch a pipeline, interactively collecting params

---

<a href="../../../../../../tools/nf_core/launch.py#L41"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Launch`

Class to hold config option to launch a pipeline

<a href="../../../../../../tools/nf_core/launch.py#L44"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline=None,
    revision=None,
    command_only=False,
    params_in=None,
    params_out=None,
    save_all=False,
    show_hidden=False,
    url=None,
    web_id=None
)
```

Initialise the Launcher class

**Args:**

- <b>`schema`</b>: An nf_core.schema.PipelineSchema() object

---

<a href="../../../../../../tools/nf_core/launch.py#L622"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_command`

```python
build_command()
```

Build the nextflow run command based on what we know

---

<a href="../../../../../../tools/nf_core/launch.py#L183"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_pipeline_schema`

```python
get_pipeline_schema()
```

Load and validate the schema from the supplied pipeline

---

<a href="../../../../../../tools/nf_core/launch.py#L309"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_web_launch_response`

```python
get_web_launch_response()
```

Given a URL for a web-gui launch response, recursively query it until results are ready.

---

<a href="../../../../../../tools/nf_core/launch.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_pipeline`

```python
launch_pipeline()
```

---

<a href="../../../../../../tools/nf_core/launch.py#L271"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_web_gui`

```python
launch_web_gui()
```

Send schema to nf-core website and launch input GUI

---

<a href="../../../../../../tools/nf_core/launch.py#L653"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_workflow`

```python
launch_workflow()
```

Launch nextflow if required

---

<a href="../../../../../../tools/nf_core/launch.py#L244"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `merge_nxf_flag_schema`

```python
merge_nxf_flag_schema()
```

Take the Nextflow flag schema and merge it with the pipeline schema

---

<a href="../../../../../../tools/nf_core/launch.py#L595"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_param_header`

```python
print_param_header(param_id, param_obj)
```

---

<a href="../../../../../../tools/nf_core/launch.py#L417"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_group`

```python
prompt_group(group_id, group_obj)
```

Prompt for edits to a group of parameters (subschema in 'definitions')

**Args:**

- <b>`group_id`</b>: Paramater ID (string)
- <b>`group_obj`</b>: JSON Schema keys (dict)

**Returns:**

- <b>`Dict of param_id`</b>: val answers

---

<a href="../../../../../../tools/nf_core/launch.py#L400"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_param`

```python
prompt_param(param_id, param_obj, is_required, answers)
```

Prompt for a single parameter

---

<a href="../../../../../../tools/nf_core/launch.py#L376"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_schema`

```python
prompt_schema()
```

Go through the pipeline schema and prompt user to change defaults

---

<a href="../../../../../../tools/nf_core/launch.py#L256"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_web_gui`

```python
prompt_web_gui()
```

Ask whether to use the web-based or cli wizard to collect params

---

<a href="../../../../../../tools/nf_core/launch.py#L350"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `sanitise_web_response`

```python
sanitise_web_response()
```

The web builder returns everything as strings. Use the functions defined in the cli wizard to convert to the correct types.

---

<a href="../../../../../../tools/nf_core/launch.py#L229"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `set_schema_inputs`

```python
set_schema_inputs()
```

Take the loaded schema and set the defaults as the input parameters If a nf_params.json file is supplied, apply these over the top

---

<a href="../../../../../../tools/nf_core/launch.py#L469"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `single_param_to_questionary`

```python
single_param_to_questionary(param_id, param_obj, answers=None, print_help=True)
```

Convert a JSONSchema param to a Questionary question

**Args:**

- <b>`param_id`</b>: Parameter ID (string)
- <b>`param_obj`</b>: JSON Schema keys (dict)
- <b>`answers`</b>: Optional preexisting answers (dict)
- <b>`print_help`</b>: If description and help_text should be printed (bool)

**Returns:**
Single Questionary dict, to be appended to questions list

---

<a href="../../../../../../tools/nf_core/launch.py#L609"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `strip_default_params`

```python
strip_default_params()
```

Strip parameters if they have not changed from the default

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
