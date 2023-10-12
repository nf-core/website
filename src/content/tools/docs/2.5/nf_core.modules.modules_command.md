<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_command`

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleCommand`

Base class for the 'nf-core modules' commands

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(dir, remote_url=None, branch=None, no_pull=False, hide_progress=False)
```

Initialise the ModulesCommand object

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L75"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `clear_module_dir`

```python
clear_module_dir(module_name, module_dir)
```

Removes all files in the module directory

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_modules`

```python
get_local_modules()
```

Get the local modules in a pipeline

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_modules_clone_modules`

```python
get_modules_clone_modules()
```

Get the modules available in a clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_modules_file`

```python
has_modules_file()
```

Checks whether a module.json file has been created and creates one if it is missing

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L55"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_valid_directory`

```python
has_valid_directory()
```

Check that we were given a pipeline or clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `install_module_files`

```python
install_module_files(module_name, module_version, modules_repo, install_dir)
```

Installs a module into the given directory

**Args:**

- <b>`module_name`</b> (str): The name of the module
- <b>`module_versioN`</b> (str): Git SHA for the version of the module to be installed
- <b>`modules_repo`</b> (ModulesRepo): A correctly configured ModulesRepo object
- <b>`install_dir`</b> (str): The path to where the module should be installed (should be the 'modules/' dir of the pipeline)

**Returns:**

- <b>`(bool)`</b>: Whether the operation was successful of not

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_lint_config`

```python
load_lint_config()
```

Parse a pipeline lint config file.

Look for a file called either `.nf-core-lint.yml` or `.nf-core-lint.yaml` in the pipeline root directory and parse it. (`.yml` takes precedence).

Add parsed config to the `self.lint_config` class attribute.

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L94"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `modules_from_repo`

```python
modules_from_repo(repo_name)
```

Gets the modules installed from a certain repository

**Args:**

- <b>`repo_name`</b> (str): The name of the repository

**Returns:**

- <b>`[str]`</b>: The names of the modules

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
