<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/components/nfcore_component.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.components.nfcore_component`

The NFCoreComponent class holds information and utility functions for a single module or subworkflow

---

<a href="../../../../../../tools/nf_core/components/nfcore_component.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `NFCoreComponent`

A class to hold the information about a nf-core module or subworkflow. Includes functionality for linting.

<a href="../../../../../../tools/nf_core/components/nfcore_component.py#L13"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    component_name,
    repo_url,
    component_dir,
    repo_type,
    base_dir,
    component_type,
    remote_component=True
)
```

Initialize the object

**Args:**

- <b>`component_name`</b> (str): The name of the module or subworkflow
- <b>`repo_url`</b> (str): The URL of the repository
- <b>`component_dir`</b> (Path): The absolute path to the module or subworkflow
- <b>`repo_type`</b> (str): Either 'pipeline' or 'modules' depending on whether the directory is a pipeline or clone of nf-core/modules.
- <b>`base_dir`</b> (Path): The absolute path to the pipeline base dir
- <b>`component_type`</b> (str): Either 'modules' or 'subworkflows'
- <b>`remote_component`</b> (bool): Whether the module is to be treated as a nf-core or local component

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
