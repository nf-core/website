<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_command`

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleCommand`

Base class for the 'nf-core modules' commands

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(dir)
```

Initialise the ModulesCommand object

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L109"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `clear_module_dir`

```python
clear_module_dir(module_name, module_dir)
```

Removes all files in the module directory

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L128"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_module_file`

```python
download_module_file(
    module_name,
    module_version,
    modules_repo,
    install_folder,
    module_dir
)
```

Downloads the files of a module from the remote repo

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L161"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `dump_modules_json`

```python
dump_modules_json(modules_json)
```

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L42"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_pipeline_modules`

```python
get_pipeline_modules()
```

Get the modules installed in the current directory.

If the current directory is a pipeline, the `module_names` field is set to a dictionary indexed by the different installation repositories in the directory. If the directory is a clone of nf-core/modules the filed is set to `{"modules": modules_in_dir}`

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L99"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_modules_file`

```python
has_modules_file()
```

Checks whether a module.json file has been created and creates one if it is missing

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L82"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_valid_directory`

```python
has_valid_directory()
```

Check that we were given a pipeline or clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L166"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_lint_config`

```python
load_lint_config()
```

Parse a pipeline lint config file.

Look for a file called either `.nf-core-lint.yml` or `.nf-core-lint.yaml` in the pipeline root directory and parse it. (`.yml` takes precedence).

Add parsed config to the `self.lint_config` class attribute.

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L143"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load_modules_json`

```python
load_modules_json()
```

Loads the modules.json file

---

<a href="../../../../../../tools/nf_core/modules/modules_command.py#L154"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update_modules_json`

```python
update_modules_json(modules_json, repo_name, module_name, module_version)
```

Updates the 'module.json' file with new module info

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
