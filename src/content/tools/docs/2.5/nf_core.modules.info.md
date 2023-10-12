<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/info.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.info`

## **Global Variables**

- **NF_CORE_MODULES_REMOTE**

---

<a href="../../../../../../tools/nf_core/modules/info.py#L23"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleInfo`

Class to print information of a module.

Attributes
---------- meta : YAML object stores the information from meta.yml file local_path : str path of the local modules remote_location : str remote repository URL local : bool indicates if the module is locally installed or not repo_type : str repository type. Can be either 'pipeline' or 'modules' modules_json : ModulesJson object contains 'modules.json' file information from a pipeline module : str name of the tool to get information from

Methods
------- init_mod_name(module) Makes sure that we have a module name get_module_info() Given the name of a module, parse meta.yml and print usage help get_local_yaml() Attempt to get the meta.yml file from a locally installed module get_remote_yaml() Attempt to get the meta.yml file from a remote repo generate_module_info_help() Take the parsed meta.yml and generate rich help

<a href="../../../../../../tools/nf_core/modules/info.py#L58"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir, tool, remote_url, branch, no_pull)
```

---

<a href="../../../../../../tools/nf_core/modules/info.py#L185"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `generate_module_info_help`

```python
generate_module_info_help()
```

Take the parsed meta.yml and generate rich help.

**Returns:**
rich renderable

---

<a href="../../../../../../tools/nf_core/modules/info.py#L129"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_yaml`

```python
get_local_yaml()
```

Attempt to get the meta.yml file from a locally installed module.

**Returns:**

- <b>`dict or bool`</b>: Parsed meta.yml found, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/info.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_info`

```python
get_module_info()
```

Given the name of a module, parse meta.yml and print usage help.

---

<a href="../../../../../../tools/nf_core/modules/info.py#L169"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_yaml`

```python
get_remote_yaml()
```

Attempt to get the meta.yml file from a remote repo.

**Returns:**

- <b>`dict or bool`</b>: Parsed meta.yml found, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/info.py#L81"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `init_mod_name`

```python
init_mod_name(module)
```

Makes sure that we have a module name before proceeding.

**Args:**

- <b>`module`</b>: str: Module name to check

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
