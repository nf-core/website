<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.module_utils`

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `path_from_remote`

```python
path_from_remote(remote_url)
```

Extracts the path from the remote URL See https://mirrors.edge.kernel.org/pub/software/scm/git/docs/git-clone.html#URLS for the possible URL patterns

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L50"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `repo_full_name_from_remote`

```python
repo_full_name_from_remote(remote_url)
```

Extracts the path from the remote URL See https://mirrors.edge.kernel.org/pub/software/scm/git/docs/git-clone.html#URLS for the possible URL patterns

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L74"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_installed_modules`

```python
get_installed_modules(dir, repo_type='modules')
```

Make a list of all modules installed in this repository

Returns a tuple of two lists, one for local modules and one for nf-core modules. The local modules are represented as direct filepaths to the module '.nf' file. Nf-core module are returned as file paths to the module directories. In case the module contains several tools, one path to each tool directory is returned.

returns (local_modules, nfcore_modules)

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_repo_type`

```python
get_repo_type(dir, repo_type=None, use_prompt=True)
```

Determine whether this is a pipeline repository or a clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L187"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_module_version_sha`

```python
prompt_module_version_sha(module, modules_repo, installed_sha=None)
```

Creates an interactive questionary prompt for selecting the module version

**Args:**

- <b>`module`</b> (str): Module name
- <b>`modules_repo`</b> (ModulesRepo): Modules repo the module originate in
- <b>`installed_sha`</b> (str): Optional extra argument to highlight the current installed version

**Returns:**

- <b>`git_sha`</b> (str): The selected version of the module

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleException`

Exception raised when there was an error with module commands

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
