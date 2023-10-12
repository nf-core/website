<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/launch.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.launch`

Launch a pipeline, interactively collecting params

---

<a href="../../../../../../tools/nf_core/launch.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `launch_pipeline`

```python
launch_pipeline(workflow, params_local_uri, direct)
```

---

<a href="../../../../../../tools/nf_core/launch.py#L44"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Launch`

Class to hold config option to launch a pipeline

<a href="../../../../../../tools/nf_core/launch.py#L47"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(workflow)
```

Initialise the class with empty placeholder vars

---

<a href="../../../../../../tools/nf_core/launch.py#L355"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_command`

```python
build_command()
```

Build the nextflow run command based on what we know

---

<a href="../../../../../../tools/nf_core/launch.py#L143"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `collect_pipeline_param_defaults`

```python
collect_pipeline_param_defaults()
```

Collect the default params and values from the workflow

---

<a href="../../../../../../tools/nf_core/launch.py#L388"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create_nfx_params_file`

```python
create_nfx_params_file()
```

---

<a href="../../../../../../tools/nf_core/launch.py#L91"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_wf`

```python
get_local_wf()
```

Check if this workflow has a local copy and use nextflow to pull it if not

---

<a href="../../../../../../tools/nf_core/launch.py#L260"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `group_parameters`

```python
group_parameters()
```

Groups parameters by their 'group' property.

**Args:**

- <b>`parameters`</b> (list): Collection of parameter objects.

**Returns:**

- <b>`dict`</b>: Parameter objects grouped by the `group` property.

---

<a href="../../../../../../tools/nf_core/launch.py#L405"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `launch_workflow`

```python
launch_workflow()
```

Launch nextflow if required

---

<a href="../../../../../../tools/nf_core/launch.py#L109"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `parse_parameter_settings`

```python
parse_parameter_settings(params_local_uri=None)
```

Load full parameter info from the pipeline parameters.settings.json file

---

<a href="../../../../../../tools/nf_core/launch.py#L222"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_core_nxf_flags`

```python
prompt_core_nxf_flags()
```

Ask the user if they want to override any default values

---

<a href="../../../../../../tools/nf_core/launch.py#L274"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_param_flags`

```python
prompt_param_flags()
```

Prompts the user for parameter input values and validates them.

---

<a href="../../../../../../tools/nf_core/launch.py#L398"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `write_params_as_full_json`

```python
write_params_as_full_json(outdir='/Users/mitochondrium/tools')
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
