<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_json`

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L23"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesJson`

An object for handling a 'modules.json' file in a pipeline

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L28"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir)
```

Initialise the object.

**Args:**

- <b>`pipeline_dir`</b> (str): The pipeline directory

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L595"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `add_patch_entry`

```python
add_patch_entry(module_name, repo_name, patch_filename, write_file=True)
```

Adds (or replaces) the patch entry for a module

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L422"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_up_to_date`

```python
check_up_to_date()
```

Checks whether the modules installed in the directory are consistent with the entries in the 'modules.json' file and vice versa.

If a module has an entry in the 'modules.json' file but is missing in the directory, we first try to reinstall the module from the remote and if that fails we remove the entry in 'modules.json'.

If a module is installed but the entry in 'modules.json' is missing we iterate through the commit log in the remote to try to determine the SHA.

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L40"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create`

```python
create()
```

Creates the modules.json file from the modules installed in the pipeline directory

**Raises:**

- <b>`UserWarning`</b>: If the creation fails

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L193"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `determine_module_branches_and_shas`

```python
determine_module_branches_and_shas(repo_name, remote_url, modules)
```

Determines what branch and commit sha each module in the pipeline belong to

Assumes all modules are installed from the default branch. If it fails to find the module in the default branch, it prompts the user with the available branches

**Args:**

- <b>`repo_name`</b> (str): The name of the module repository
- <b>`remote_url`</b> (str): The url to the remote repository
- <b>`modules`</b> ([str]): List of names of installed modules from the repository

**Returns:**

- <b>`(dict[str, dict[str, str]])`</b>: The module.json entries for the modules from the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `dir_tree_uncovered`

```python
dir_tree_uncovered(modules_dir, repos)
```

Does a BFS of the modules directory to look for directories that are not tracked by a remote. The 'repos' argument contains the directories that are currently covered by remote, and it and its subdirectories are therefore ignore.

**Args:**

- <b>`module_dir`</b> (Path): Base path of modules in pipeline
- <b>`repos`</b> ([ Path ]): List of repos that are covered by a remote

**Returns:**

- <b>`dirs_not_covered`</b> ([ Path ]): A list of directories that are currently not covered by any remote.

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L769"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `dump`

```python
dump()
```

Sort the modules.json, and write it to file

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L275"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_correct_commit_sha`

```python
find_correct_commit_sha(module_name, module_path, modules_repo)
```

Returns the SHA for the latest commit where the local files are identical to the remote files

**Args:**

- <b>`module_name`</b> (str): Name of module
- <b>`module_path`</b> (str): Path to module in local repo
- <b>`module_repo`</b> (str): Remote repo for module

**Returns:**

- <b>`commit_sha`</b> (str): The latest commit SHA where local files are identical to remote files, or None if no commit is found

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L732"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_all_modules`

```python
get_all_modules()
```

Retrieves all pipeline modules that are reported in the modules.json

**Returns:**

- <b>`(dict[str, [str]])`</b>: Dictionary indexed with the repo names, with a list of modules as values

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L718"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_git_url`

```python
get_git_url(repo_name)
```

Returns the git url of a repo

**Args:**

- <b>`repo_name`</b> (str): Name of the repository

**Returns:**

- <b>`(str)`</b>: The git url of the repository if it exists, None otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L750"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_branch`

```python
get_module_branch(module, repo_name)
```

Gets the branch from which the module was installed

**Returns:**

- <b>`(str)`</b>: The branch name

**Raises:**

- <b>`LookupError`</b>: If their is no branch entry in the `modules.json`

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L697"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_version`

```python
get_module_version(module_name, repo_name)
```

Returns the version of a module

**Args:**

- <b>`module_name`</b> (str): Name of the module
- <b>`repo_name`</b> (str): Name of the repository

**Returns:**

- <b>`(str)`</b>: The git SHA of the module if it exists, None otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L686"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_modules_json`

```python
get_modules_json()
```

Returns a copy of the loaded modules.json

**Returns:**

- <b>`(dict)`</b>: A copy of the loaded modules.json

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L609"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_patch_fn`

```python
get_patch_fn(module_name, repo_name)
```

Get the patch filename of a module

**Args:**

- <b>`module_name`</b> (str): The name of the module
- <b>`repo_name`</b> (str): The name of the repository containing the module

**Returns:**

- <b>`(str)`</b>: The patch filename for the module, None if not present

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L86"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_pipeline_module_repositories`

```python
get_pipeline_module_repositories(modules_dir, repos=None)
```

Finds all module repositories in the modules directory. Ignores the local modules.

**Args:**

- <b>`modules_dir`</b> (Path): base directory for the module files Returns repos ([ (str, str, str) ]),
- <b>`renamed_dirs`</b> (dict[Path, Path]): List of tuples of repo name, repo remote URL and path to modules in repo

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L364"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_git_url_and_modules`

```python
has_git_url_and_modules()
```

Check that all repo entries in the modules.json has a git url and a modules dict entry

**Returns:**

- <b>`(bool)`</b>: True if they are found for all repos, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L523"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load`

```python
load()
```

Loads the modules.json file into the variable 'modules_json'

Sets the modules_json attribute to the loaded file.

**Raises:**

- <b>`UserWarning`</b>: If the modules.json file is not found

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L673"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `module_present`

```python
module_present(module_name, repo_name)
```

Checks if a module is present in the modules.json file

**Args:**

- <b>`module_name`</b> (str): Name of the module
- <b>`repo_name`</b> (str): Name of the repository

**Returns:**

- <b>`(bool)`</b>: Whether the module is present in the 'modules.json' file

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L295"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `move_module_to_local`

```python
move_module_to_local(module, repo_name)
```

Move a module to the 'local' directory

**Args:**

- <b>`module`</b> (str): The name of the modules
- <b>`repo_name`</b> (str): The name of the repository the module resides in

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L385"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `reinstall_repo`

```python
reinstall_repo(repo_name, remote_url, module_entries)
```

Reinstall modules from a repository

**Args:**

- <b>`repo_name`</b> (str): The name of the repository
- <b>`remote_url`</b> (str): The git url of the remote repository
- <b>`modules`</b> ([ dict[str, dict[str, str]] ]): Module entries with branch and git sha info

**Returns:**

- <b>`([ str ])`</b>: List of modules that we failed to install

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L567"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_entry`

```python
remove_entry(module_name, repo_name)
```

Removes an entry from the 'modules.json' file.

**Args:**

- <b>`module_name`</b> (str): Name of the module to be removed
- <b>`repo_name`</b> (str): Name of the repository containing the module

**Returns:**

- <b>`(bool)`</b>: True if the removal was successful, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L661"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `repo_present`

```python
repo_present(repo_name)
```

Checks if a repo is present in the modules.json file

**Args:**

- <b>`repo_name`</b> (str): Name of the repository

**Returns:**

- <b>`(bool)`</b>: Whether the repo exists in the modules.json

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L625"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `try_apply_patch_reverse`

```python
try_apply_patch_reverse(module, repo_name, patch_relpath, module_dir)
```

Try reverse applying a patch file to the modified module files

**Args:**

- <b>`module`</b> (str): The name of the module
- <b>`repo_name`</b> (str): The name of the repository where the module resides
- <b>`patch_relpath`</b> (Path | str): The path to patch file in the pipeline
- <b>`module_dir`</b> (Path | str): The module directory in the pipeline

**Returns:**

- <b>`(Path | str)`</b>: The path of the folder where the module patched files are

**Raises:**

- <b>`LookupError`</b>: If patch was not applied

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L316"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `unsynced_modules`

```python
unsynced_modules()
```

Compute the difference between the modules in the directory and the modules in the 'modules.json' file. This is done by looking at all directories containing a 'main.nf' file

**Returns:**

- <b>`(untrack_dirs ([ Path ]), missing_installation (dict))`</b>: Directories that are not tracked by the modules.json file, and modules in the modules.json where the installation directory is missing

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L539"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update`

```python
update(modules_repo, module_name, module_version, write_file=True)
```

Updates the 'module.json' file with new module info

**Args:**

- <b>`modules_repo`</b> (ModulesRepo): A ModulesRepo object configured for the new module
- <b>`module_name`</b> (str): Name of new module
- <b>`module_version`</b> (str): git SHA for the new module entry
- <b>`write_file`</b> (bool): whether to write the updated modules.json to a file.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
