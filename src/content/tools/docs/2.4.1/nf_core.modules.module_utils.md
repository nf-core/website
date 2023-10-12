<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.module_utils`

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `module_exist_in_repo`

```python
module_exist_in_repo(module_name, modules_repo)
```

Checks whether a module exists in a branch of a GitHub repository

**Args:**

- <b>`module_name`</b> (str): Name of module
- <b>`modules_repo`</b> (ModulesRepo): A ModulesRepo object configured for the repository in question

**Returns:**

- <b>`boolean`</b>: Whether the module exist in the repo or not.

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L43"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_module_git_log`

```python
get_module_git_log(
    module_name,
    modules_repo=None,
    per_page=30,
    page_nbr=1,
    since='2021-07-07T00:00:00Z'
)
```

Fetches the commit history the of requested module since a given date. The default value is not arbitrary - it is the last time the structure of the nf-core/modules repository was had an update breaking backwards compatibility.

**Args:**

- <b>`module_name`</b> (str): Name of module
- <b>`modules_repo`</b> (ModulesRepo): A ModulesRepo object configured for the repository in question
- <b>`per_page`</b> (int): Number of commits per page returned by API
- <b>`page_nbr`</b> (int): Page number of the retrieved commits
- <b>`since`</b> (str): Only show commits later than this timestamp.
- <b>`Time should be given in ISO-8601 format`</b>: YYYY-MM-DDTHH:MM:SSZ.

**Returns:**

- <b>`[ dict ]`</b>: List of commit SHAs and associated (truncated) message

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L90"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_commit_info`

```python
get_commit_info(commit_sha, repo_name='nf-core/modules')
```

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L124"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `create_modules_json`

```python
create_modules_json(pipeline_dir)
```

Create the modules.json files

**Args:**

- <b>`pipeline_dir`</b> (str): The directory where the `modules.json` should be created

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L197"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `find_correct_commit_sha`

```python
find_correct_commit_sha(module_name, module_path, modules_repo)
```

Returns the SHA for the latest commit where the local files are identical to the remote files

**Args:**

- <b>`module_name`</b> (str): Name of module
- <b>`module_path`</b> (str): Path to module in local repo
- <b>`module_repo`</b> (str): Remote repo for module

**Returns:**

- <b>`commit_sha`</b> (str): The latest commit SHA where local files are identical to remote files

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L225"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `iterate_commit_log_page`

```python
iterate_commit_log_page(module_name, module_path, modules_repo, commit_shas)
```

Iterates through a list of commits for a module and checks if the local file contents match the remote

**Args:**

- <b>`module_name`</b> (str): Name of module
- <b>`module_path`</b> (str): Path to module in local repo
- <b>`module_repo`</b> (str): Remote repo for module
- <b>`commit_shas`</b> ([ str ]): List of commit SHAs for module, sorted in descending order

**Returns:**

- <b>`commit_sha`</b> (str): The latest commit SHA from 'commit_shas' where local files are identical to remote files

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L252"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `local_module_equal_to_commit`

```python
local_module_equal_to_commit(local_files, module_name, modules_repo, commit_sha)
```

Compares the local module files to the module files for the given commit sha

**Args:**

- <b>`local_files`</b> ([ str ]): Contents of local files. `None` if files doesn't exist
- <b>`module_name`</b> (str): Name of module
- <b>`module_repo`</b> (str): Remote repo for module
- <b>`commit_sha`</b> (str): Commit SHA for remote version to compare against local version

**Returns:**

- <b>`bool`</b>: Whether all local files are identical to remote version

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L293"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_installed_modules`

```python
get_installed_modules(dir, repo_type='modules')
```

Make a list of all modules installed in this repository

Returns a tuple of two lists, one for local modules and one for nf-core modules. The local modules are represented as direct filepaths to the module '.nf' file. Nf-core module are returned as file paths to the module directories. In case the module contains several tools, one path to each tool directory is returned.

returns (local_modules, nfcore_modules)

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L349"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_repo_type`

```python
get_repo_type(dir, repo_type=None, use_prompt=True)
```

Determine whether this is a pipeline repository or a clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L409"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `verify_pipeline_dir`

```python
verify_pipeline_dir(dir)
```

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L439"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_module_version_sha`

```python
prompt_module_version_sha(module, modules_repo, installed_sha=None)
```

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L483"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `sha_exists`

```python
sha_exists(sha, modules_repo)
```

---

<a href="../../../../../../tools/nf_core/modules/module_utils.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleException`

Exception raised when there was an error with module commands

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
