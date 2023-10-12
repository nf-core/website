<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/nfcore_module.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.nfcore_module`

The NFCoreModule class holds information and utility functions for a single module

---

<a href="../../../../../../tools/nf_core/modules/nfcore_module.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `NFCoreModule`

A class to hold the information about a nf-core module Includes functionality for linting

<a href="../../../../../../tools/nf_core/modules/nfcore_module.py#L13"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    module_name,
    repo_url,
    module_dir,
    repo_type,
    base_dir,
    nf_core_module=True
)
```

Initialize the object

**Args:**

- <b>`module_dir`</b> (Path): The absolute path to the module
- <b>`repo_type`</b> (str): Either 'pipeline' or 'modules' depending on whether the directory is a pipeline or clone of nf-core/modules.
- <b>`base_dir`</b> (Path): The absolute path to the pipeline base dir
- <b>`nf_core_module`</b> (bool): Whether the module is to be treated as a nf-core or local module

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
