<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/workflow/workflow.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.workflow.workflow`

---

<a href="../../../../../../tools/nf_core/workflow/workflow.py#L6"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Workflow`

nf-core workflow object that holds run parameter information.

**Args:**

- <b>`name`</b> (str): Workflow name.
- <b>`parameters_json`</b> (str): Workflow parameter data in JSON.

<a href="../../../../../../tools/nf_core/workflow/workflow.py#L13"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(name, parameters_json)
```

---

<a href="../../../../../../tools/nf_core/workflow/workflow.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `in_full_json`

```python
in_full_json(indent=0)
```

Converts the Parameter list in a complete parameter JSON for schema validation.

**Returns:**

- <b>`str`</b>: JSON formatted parameters.

---

<a href="../../../../../../tools/nf_core/workflow/workflow.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `in_nextflow_json`

```python
in_nextflow_json(indent=0)
```

Converts the Parameter list in a workflow readable parameter JSON file.

**Returns:**

- <b>`str`</b>: JSON formatted parameters.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
