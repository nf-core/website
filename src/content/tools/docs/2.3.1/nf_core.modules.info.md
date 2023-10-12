<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/info.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.info`

---

<a href="../../../../../../tools/nf_core/modules/info.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleInfo`

<a href="../../../../../../tools/nf_core/modules/info.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir, tool)
```

---

<a href="../../../../../../tools/nf_core/modules/info.py#L115"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `generate_module_info_help`

```python
generate_module_info_help()
```

Take the parsed meta.yml and generate rich help.

**Returns:**
rich renderable

---

<a href="../../../../../../tools/nf_core/modules/info.py#L57"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_yaml`

```python
get_local_yaml()
```

Attempt to get the meta.yml file from a locally installed module.

**Returns:**

- <b>`dict or bool`</b>: Parsed meta.yml found, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/info.py#L40"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_info`

```python
get_module_info()
```

Given the name of a module, parse meta.yml and print usage help.

---

<a href="../../../../../../tools/nf_core/modules/info.py#L85"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_yaml`

```python
get_remote_yaml()
```

Attempt to get the meta.yml file from a remote repo.

**Returns:**

- <b>`dict or bool`</b>: Parsed meta.yml found, False otherwise

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
