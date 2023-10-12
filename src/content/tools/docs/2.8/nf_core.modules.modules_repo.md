<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_repo`

## **Global Variables**

- **NFCORE_DIR**
- **NF_CORE_MODULES_NAME**
- **NF_CORE_MODULES_REMOTE**
- **NF_CORE_MODULES_DEFAULT_BRANCH**

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `RemoteProgressbar`

An object to create a progressbar for when doing an operation with the remote. Note that an initialized rich Progress (progress bar) object must be past during initialization.

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L31"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(progress_bar, repo_name, remote_url, operation)
```

Initializes the object and adds a task to the progressbar passed as 'progress_bar'

**Args:**

- <b>`progress_bar`</b> (rich.progress.Progress): A rich progress bar object
- <b>`repo_name`</b> (str): Name of the repository the operation is performed on
- <b>`remote_url`</b> (str): Git URL of the repository the operation is performed on
- <b>`operation`</b> (str): The operation performed on the repository, i.e. 'Pulling', 'Cloning' etc.

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L49"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update`

```python
update(op_code, cur_count, max_count=None, message='')
```

Overrides git.RemoteProgress.update. Called every time there is a change in the remote operation

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L61"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesRepo`

An object to store details about the repository being used for modules.

Used by the `nf-core modules` top-level command with -r and -b flags, so that this can be used in the same way by all sub-commands.

We keep track of the pull-status of the different installed repos in the static variable local_repo_status. This is so we don't need to pull a remote several times in one command.

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L115"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(remote_url=None, branch=None, no_pull=False, hide_progress=False)
```

Initializes the object and clones the git repository if it is not already present

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L268"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `branch_exists`

```python
branch_exists()
```

Verifies that the branch exists in the repository by trying to check it out

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L296"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `checkout`

```python
checkout(commit)
```

Checks out the repository at the requested commit

**Args:**

- <b>`commit`</b> (str): Git SHA of the commit

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L290"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `checkout_branch`

```python
checkout_branch()
```

Checks out the specified branch of the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L305"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `component_exists`

```python
component_exists(component_name, component_type, checkout=True, commit=None)
```

Check if a module/subworkflow exists in the branch of the repo

**Args:**

- <b>`component_name`</b> (str): The name of the module/subworkflow

**Returns:**

- <b>`(bool)`</b>: Whether the module/subworkflow exists in this branch of the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L452"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_avail_components`

```python
get_avail_components(component_type, checkout=True, commit=None)
```

Gets the names of the modules/subworkflows in the repository. They are detected by checking which directories have a 'main.nf' file

**Returns:**

- <b>`([ str ])`</b>: The module/subworkflow names

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L432"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_commit_info`

```python
get_commit_info(sha)
```

Fetches metadata about the commit (dates, message, etc.)

**Args:**

- <b>`commit_sha`</b> (str): The SHA of the requested commit

**Returns:**

- <b>`message`</b> (str): The commit message for the requested commit
- <b>`date`</b> (str): The commit date for the requested commit

**Raises:**

- <b>`LookupError`</b>: If the search for the commit fails

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L317"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_component_dir`

```python
get_component_dir(component_name, component_type)
```

Returns the file path of a module/subworkflow directory in the repo. Does not verify that the path exists.

**Args:**

- <b>`component_name`</b> (str): The name of the module/subworkflow

**Returns:**

- <b>`component_path`</b> (str): The path of the module/subworkflow in the local copy of the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L390"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_component_git_log`

```python
get_component_git_log(component_name, component_type, depth=None)
```

Fetches the commit history the of requested module/subworkflow since a given date. The default value is not arbitrary - it is the last time the structure of the nf-core/modules repository was had an update breaking backwards compatibility.

**Args:**

- <b>`component_name`</b> (str): Name of module/subworkflow
- <b>`modules_repo`</b> (ModulesRepo): A ModulesRepo object configured for the repository in question

**Returns:**

- <b>`( dict )`</b>: Iterator of commit SHAs and associated (truncated) message

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L260"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_default_branch`

```python
get_default_branch()
```

Gets the default branch for the repo (the branch origin/HEAD is pointing to)

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L419"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_latest_component_version`

```python
get_latest_component_version(component_name, component_type)
```

Returns the latest commit in the repository

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L477"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_meta_yml`

```python
get_meta_yml(component_type, module_name)
```

Returns the contents of the 'meta.yml' file of a module

**Args:**

- <b>`module_name`</b> (str): The name of the module

**Returns:**

- <b>`(str)`</b>: The contents of the file in text format

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L90"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_branches`

```python
get_remote_branches(remote_url)
```

Get all branches from a remote repository

**Args:**

- <b>`remote_url`</b> (str): The git url to the remote repository

**Returns:**

- <b>`(set[str])`</b>: All branches found in the remote

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L332"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `install_component`

```python
install_component(component_name, install_dir, commit, component_type)
```

Install the module/subworkflow files into a pipeline at the given commit

**Args:**

- <b>`component_name`</b> (str): The name of the module/subworkflow
- <b>`install_dir`</b> (str): The path where the module/subworkflow should be installed
- <b>`commit`</b> (str): The git SHA for the version of the module/subworkflow to be installed

**Returns:**

- <b>`(bool)`</b>: Whether the operation was successful or not

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L76"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `local_repo_synced`

```python
local_repo_synced(repo_name)
```

Checks whether a local repo has been cloned/pull in the current session

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L364"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `module_files_identical`

```python
module_files_identical(module_name, base_path, commit)
```

Checks whether the module files in a pipeline are identical to the ones in the remote

**Args:**

- <b>`module_name`</b> (str): The name of the module
- <b>`base_path`</b> (str): The path to the module in the pipeline

**Returns:**

- <b>`(bool)`</b>: Whether the pipeline files are identical to the repo files

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L240"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_branch`

```python
setup_branch(branch)
```

Verify that we have a branch and otherwise use the default one. The branch is then checked out to verify that it exists in the repo.

**Args:**

- <b>`branch`</b> (str): Name of branch

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L169"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_local_repo`

```python
setup_local_repo(remote, branch, hide_progress=True)
```

Sets up the local git repository. If the repository has been cloned previously, it returns a git.Repo object of that clone. Otherwise it tries to clone the repository from the provided remote URL and returns a git.Repo of the new clone.

**Args:**

- <b>`remote`</b> (str): git url of remote
- <b>`branch`</b> (str): name of branch to use Sets self.repo

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L425"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `sha_exists_on_branch`

```python
sha_exists_on_branch(sha)
```

Verifies that a given commit sha exists on the branch

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L83"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update_local_repo_status`

```python
update_local_repo_status(repo_name, up_to_date)
```

Updates the clone/pull status of a local repo

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L277"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `verify_branch`

```python
verify_branch()
```

Verifies the active branch conforms do the correct directory structure

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L149"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `verify_sha`

```python
verify_sha(prompt, sha)
```

Verify that 'sha' and 'prompt' arguments are not provided together. Verify that the provided SHA exists in the repo.

**Arguments:**

- <b>`prompt`</b> (bool): prompt asking for SHA
- <b>`sha`</b> (str): provided sha

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
