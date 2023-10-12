<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_json`

## **Global Variables**

- **NF_CORE_MODULES_NAME**
- **NF_CORE_MODULES_REMOTE**

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L28"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesJson`

An object for handling a 'modules.json' file in a pipeline

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir)
```

Initialise the object.

**Args:**

- <b>`pipeline_dir`</b> (str): The pipeline directory

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L738"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `add_patch_entry`

```python
add_patch_entry(
    module_name,
    repo_url,
    install_dir,
    patch_filename,
    write_file=True
)
```

Adds (or replaces) the patch entry for a module

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L533"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_up_to_date`

```python
check_up_to_date()
```

Checks whether the modules and subworkflows installed in the directory are consistent with the entries in the 'modules.json' file and vice versa.

If a module/subworkflow has an entry in the 'modules.json' file but is missing in the directory, we first try to reinstall the module/subworkflow from the remote and if that fails we remove the entry in 'modules.json'.

If a module/subworkflow is installed but the entry in 'modules.json' is missing we iterate through the commit log in the remote to try to determine the SHA.

Check that we have the "installed_by" value in 'modules.json', otherwise add it. Assume that the modules/subworkflows were installed by an nf-core command (don't track installed by subworkflows).

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L57"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create`

```python
create()
```

Creates the modules.json file from the modules and subworkflows installed in the pipeline directory

**Raises:**

- <b>`UserWarning`</b>: If the creation fails

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L239"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `determine_branches_and_shas`

```python
determine_branches_and_shas(component_type, install_dir, remote_url, components)
```

Determines what branch and commit sha each module/subworkflow in the pipeline belongs to

Assumes all modules/subworkflows are installed from the default branch. If it fails to find the module/subworkflow in the default branch, it prompts the user with the available branches

**Args:**

- <b>`install_dir`</b> (str): The name of the directory inside modules or subworkflows where components are installed
- <b>`remote_url`</b> (str): The url to the remote repository
- <b>`components`</b> ([str]): List of names of installed modules/subworkflows from the repository

**Returns:**

- <b>`(dict[str, dict[str, str]])`</b>: The module.json entries for the modules/subworkflows from the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L204"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `dir_tree_uncovered`

```python
dir_tree_uncovered(components_directory, repos)
```

Does a BFS of the modules/subworkflos directory to look for directories that are not tracked by a remote. The 'repos' argument contains the directories that are currently covered by remote, and it and its subdirectories are therefore ignore.

**Args:**

- <b>`components_directory`</b> (Path): Base path of modules or subworkflows in pipeline
- <b>`repos`</b> ([ Path ]): List of repos that are covered by a remote

**Returns:**

- <b>`dirs_not_covered`</b> ([ Path ]): A list of directories that are currently not covered by any remote.

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L1038"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `dump`

```python
dump()
```

Sort the modules.json, and write it to file

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L337"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_correct_commit_sha`

```python
find_correct_commit_sha(
    component_type,
    component_name,
    component_path,
    modules_repo
)
```

Returns the SHA for the latest commit where the local files are identical to the remote files

**Args:**

- <b>`component_type`</b> (str): modules or subworkflows
- <b>`component_name`</b> (str): Name of module/subowrkflow
- <b>`component_path`</b> (str): Path to module/subworkflow in local repo
- <b>`modules_repo`</b> (str): Remote repo for module/subworkflow

**Returns:**

- <b>`commit_sha`</b> (str): The latest commit SHA where local files are identical to remote files, or None if no commit is found

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L930"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_all_components`

```python
get_all_components(component_type)
```

Retrieves all pipeline modules/subworkflows that are reported in the modules.json

**Returns:**

- <b>`(dict[str, [(str, str)]])`</b>: Dictionary indexed with the repo urls, with a list of tuples (component_dir, components) as values

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L1012"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_component_branch`

```python
get_component_branch(component_type, component, repo_url, install_dir)
```

Gets the branch from which the module/subworkflow was installed

**Returns:**

- <b>`(str)`</b>: The branch name

**Raises:**

- <b>`LookupError`</b>: If there is no branch entry in the `modules.json`

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L100"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_component_names_from_repo`

```python
get_component_names_from_repo(repos, directory)
```

Get component names from repositories in a pipeline.

**Args:**

- <b>`repos`</b> (list): list of repository urls
- <b>`directory`</b> (str): modules directory or subworkflows directory

**Returns:**

- <b>`[(str),[(str),(str)]]`</b>: list of tuples with repository url, component names and install directory

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L861"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_component_version`

```python
get_component_version(component_type, component_name, repo_url, install_dir)
```

Returns the version of a module or subworkflow

**Args:**

- <b>`component_name`</b> (str): Name of the module/subworkflow
- <b>`repo_url`</b> (str): URL of the repository
- <b>`install_dir`</b> (str): Name of the directory where modules/subworkflows are installed

**Returns:**

- <b>`(str)`</b>: The git SHA of the module/subworkflow if it exists, None otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L949"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_dependent_components`

```python
get_dependent_components(
    component_type,
    name,
    repo_url,
    install_dir,
    dependent_components
)
```

Retrieves all pipeline modules/subworkflows that are reported in the modules.json as being installed by the given component

**Args:**

- <b>`component_type`</b> (str): Type of component [modules, subworkflows]
- <b>`name`</b> (str): Name of the component to find dependencies for
- <b>`repo_url`</b> (str): URL of the repository containing the components
- <b>`install_dir`</b> (str): Name of the directory where components are installed

**Returns:**

- <b>`(dict[str`</b>: str,]): Dictionary indexed with the component names, with component_type as value

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L988"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_installed_by_entries`

```python
get_installed_by_entries(component_type, name)
```

Retrieves all entries of installed_by for a given component

**Args:**

- <b>`component_type`</b> (str): Type of component [modules, subworkflows]
- <b>`name`</b> (str): Name of the component to find dependencies for

**Returns:**

- <b>`(list)`</b>: The list of installed_by entries

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L884"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_version`

```python
get_module_version(module_name, repo_url, install_dir)
```

Returns the version of a module

**Args:**

- <b>`module_name`</b> (str): Name of the module
- <b>`repo_url`</b> (str): URL of the repository
- <b>`install_dir`</b> (str): Name of the directory where modules are installed

**Returns:**

- <b>`(str)`</b>: The git SHA of the module if it exists, None otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L850"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_modules_json`

```python
get_modules_json()
```

Returns a copy of the loaded modules.json

**Returns:**

- <b>`(dict)`</b>: A copy of the loaded modules.json

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L762"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_patch_fn`

```python
get_patch_fn(module_name, repo_url, install_dir)
```

Get the patch filename of a module

**Args:**

- <b>`module_name`</b> (str): The name of the module
- <b>`repo_url`</b> (str): The URL of the repository containing the module
- <b>`install_dir`</b> (str): The name of the directory where modules are installed

**Returns:**

- <b>`(str)`</b>: The patch filename for the module, None if not present

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L126"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_pipeline_module_repositories`

```python
get_pipeline_module_repositories(component_type, directory, repos=None)
```

Finds all module repositories in the modules and subworkflows directory. Ignores the local modules/subworkflows.

**Args:**

- <b>`component_type`</b> (str): modules or subworkflows
- <b>`directory`</b> (Path): base directory for the module files Returns repos ([ (str, str, str) ]),
- <b>`renamed_dirs`</b> (dict[Path, Path]): List of tuples of repo name, repo remote URL and path to modules in repo

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L907"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_subworkflow_version`

```python
get_subworkflow_version(subworkflow_name, repo_url, install_dir)
```

Returns the version of a subworkflow

**Args:**

- <b>`subworkflow_name`</b> (str): Name of the module
- <b>`repo_url`</b> (str): URL of the repository
- <b>`install_dir`</b> (str): Name of the directory where subworkflows are installed

**Returns:**

- <b>`(str)`</b>: The git SHA of the subworkflow if it exists, None otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L465"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `has_git_url_and_modules`

```python
has_git_url_and_modules()
```

Check that all repo entries in the modules.json has a git url and a modules dict entry

**Returns:**

- <b>`(bool)`</b>: True if they are found for all repos, False otherwise

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L614"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `load`

```python
load()
```

Loads the modules.json file into the variable 'modules_json'

Sets the modules_json attribute to the loaded file.

**Raises:**

- <b>`UserWarning`</b>: If the modules.json file is not found

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L834"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `module_present`

```python
module_present(module_name, repo_url, install_dir)
```

Checks if a module is present in the modules.json file

**Args:**

- <b>`module_name`</b> (str): Name of the module
- <b>`repo_url`</b> (str): URL of the repository
- <b>`install_dir`</b> (str): Name of the directory where modules are installed

**Returns:**

- <b>`(bool)`</b>: Whether the module is present in the 'modules.json' file

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L365"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `move_component_to_local`

```python
move_component_to_local(component_type, component, repo_name)
```

Move a module/subworkflow to the 'local' directory

**Args:**

- <b>`component`</b> (str): The name of the module/subworkflow
- <b>`repo_name`</b> (str): The name of the repository the module resides in

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L424"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `parse_dirs`

```python
parse_dirs(dirs, missing_installation, component_type)
```

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L1158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `recreate_dependencies`

```python
recreate_dependencies(repo, org, subworkflow)
```

Try to recreate the installed_by entries for subworkflows. Remove self installation entry from dependencies, assuming that the modules.json has been freshly created, i.e., no module or subworkflow has been installed by the user in the meantime

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L494"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `reinstall_repo`

```python
reinstall_repo(install_dir, remote_url, module_entries)
```

Reinstall modules from a repository

**Args:**

- <b>`install_dir`</b> (str): The name of directory where modules are installed
- <b>`remote_url`</b> (str): The git url of the remote repository
- <b>`module_entries`</b> ([ dict[str, dict[str, str]] ]): Module entries with branch and git sha info

**Returns:**

- <b>`([ str ])`</b>: List of modules that we failed to install

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L689"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_entry`

```python
remove_entry(component_type, name, repo_url, install_dir, removed_by=None)
```

Removes an entry from the 'modules.json' file.

**Args:**

- <b>`component_type`</b> (str): Type of component [modules, subworkflows]
- <b>`name`</b> (str): Name of the component to be removed
- <b>`repo_url`</b> (str): URL of the repository containing the component
- <b>`install_dir`</b> (str): Name of the directory where components are installed
- <b>`removed_by`</b> (str): Name of the component that wants to remove the component

**Returns:**

- <b>`(bool)`</b>: return True if the component was removed, False if it was not found or is still depended on

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L752"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_patch_entry`

```python
remove_patch_entry(module_name, repo_url, install_dir, write_file=True)
```

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L822"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L1084"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `resolve_missing_from_modules_json`

```python
resolve_missing_from_modules_json(missing_from_modules_json, component_type)
```

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L1046"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `resolve_missing_installation`

```python
resolve_missing_installation(missing_installation, component_type)
```

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L786"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L390"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `unsynced_components`

```python
unsynced_components()
```

Compute the difference between the modules/subworkflows in the directory and the modules/subworkflows in the 'modules.json' file. This is done by looking at all directories containing a 'main.nf' file

**Returns:**

- <b>`(untrack_dirs ([ Path ]), missing_installation (dict))`</b>: Directories that are not tracked by the modules.json file, and modules/subworkflows in the modules.json where the installation directory is missing

---

<a href="../../../../../../tools/nf_core/modules/modules_json.py#L633"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update`

```python
update(
    component_type,
    modules_repo,
    component_name,
    component_version,
    installed_by,
    installed_by_log=None,
    write_file=True
)
```

Updates the 'module.json' file with new module/subworkflow info

**Args:**

- <b>`component_type`</b> (str): modules or subworkflows
- <b>`modules_repo`</b> (ModulesRepo): A ModulesRepo object configured for the new module/subworkflow
- <b>`component_name`</b> (str): Name of new module/subworkflow
- <b>`component_version`</b> (str): git SHA for the new module/subworkflow entry
- <b>`installed_by_log`</b> (list): previous tracing of installed_by that needs to be added to 'modules.json'
- <b>`write_file`</b> (bool): whether to write the updated modules.json to a file.

**Returns:**

- <b>`bool`</b>: True if the module/subworkflow was successfully added to the 'modules.json' file

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
