<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/update.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.update`

---

<a href="../../../../../../tools/nf_core/modules/update.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleUpdate`

<a href="../../../../../../tools/nf_core/modules/update.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline_dir,
    force=False,
    prompt=False,
    sha=None,
    update_all=False,
    show_diff=None,
    save_diff_fn=None,
    remote_url=None,
    branch=None,
    no_pull=False
)
```

---

<a href="../../../../../../tools/nf_core/modules/update.py#L323"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_all_modules_info`

```python
get_all_modules_info(branch=None)
```

Collects the module repository, version and sha for all modules.

Information about the module version in the '.nf-core.yml' overrides the '--sha' option.

**Returns:**

- <b>`[(ModulesRepo, str, str)]`</b>: A list of tuples containing a ModulesRepo object, the module name, and the module version.

---

<a href="../../../../../../tools/nf_core/modules/update.py#L246"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_single_module_info`

```python
get_single_module_info(module)
```

Collects the module repository, version and sha for a module.

Information about the module version in the '.nf-core.yml' overrides the '--sha' option

**Args:**

- <b>`module_name`</b> (str): The name of the module to get info for.

**Returns:**

- <b>`(ModulesRepo, str, str)`</b>: The modules repo containing the module, the module name, and the module version.

**Raises:**

- <b>`LookupError`</b>: If the module is not found either in the pipeline or the modules repo.
- <b>`UserWarning`</b>: If the '.nf-core.yml' entry is not valid.

---

<a href="../../../../../../tools/nf_core/modules/update.py#L494"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `move_files_from_tmp_dir`

```python
move_files_from_tmp_dir(
    module,
    module_dir,
    install_folder,
    repo_name,
    new_version
)
```

Move the files from the temporary to the installation directory.

**Args:**

- <b>`module`</b> (str): The module name.
- <b>`module_dir`</b> (str): The path to the module directory.
- <b>`install_folder [str]`</b>: The path to the temporary installation directory.
- <b>`modules_repo`</b> (ModulesRepo): The ModulesRepo object from which the module was installed.
- <b>`new_version`</b> (str): The version of the module that was installed.

---

<a href="../../../../../../tools/nf_core/modules/update.py#L465"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_diff_file`

```python
setup_diff_file()
```

Sets up the diff file.

If the save diff option was chosen interactively, the user is asked to supply a name for the diff file.

Then creates the file for saving the diff.

---

<a href="../../../../../../tools/nf_core/modules/update.py#L519"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `try_apply_patch`

```python
try_apply_patch(
    module,
    repo_name,
    patch_relpath,
    module_dir,
    module_install_dir
)
```

Try applying a patch file to the new module files

**Args:**

- <b>`module`</b> (str): The name of the module
- <b>`repo_name`</b> (str): The name of the repository where the module resides
- <b>`patch_relpath`</b> (Path | str): The path to patch file in the pipeline
- <b>`module_dir`</b> (Path | str): The module directory in the pipeline
- <b>`module_install_dir`</b> (Path | str): The directory where the new module file have been installed

**Returns:**

- <b>`(bool)`</b>: Whether the patch application was successful

---

<a href="../../../../../../tools/nf_core/modules/update.py#L69"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update`

```python
update(module=None)
```

Updates a specified module or all modules modules in a pipeline.

**Args:**

- <b>`module`</b> (str): The name of the module to update.

**Returns:**

- <b>`bool`</b>: True if the update was successful, False otherwise.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
