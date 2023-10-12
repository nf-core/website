<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/sync.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.sync`

Synchronise a pipeline TEMPLATE branch with the template.

---

<a href="../../../../../../tools/nf_core/sync.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SyncException`

Exception raised when there was an error with TEMPLATE branch synchronisation

---

<a href="../../../../../../tools/nf_core/sync.py#L30"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PullRequestException`

Exception raised when there was an error creating a Pull-Request on GitHub.com

---

<a href="../../../../../../tools/nf_core/sync.py#L36"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineSync`

Object to hold syncing information and results.

**Args:**

- <b>`pipeline_dir`</b> (str): The path to the Nextflow pipeline root directory
- <b>`from_branch`</b> (str): The branch to use to fetch config vars. If not set, will use current active branch
- <b>`make_pr`</b> (bool): Set this to `True` to create a GitHub pull-request with the changes
- <b>`gh_username`</b> (str): GitHub username
- <b>`gh_repo`</b> (str): GitHub repository name

**Attributes:**

- <b>`pipeline_dir`</b> (str): Path to target pipeline directory
- <b>`from_branch`</b> (str): Repo branch to use when collecting workflow variables. Default: active branch.
- <b>`original_branch`</b> (str): Repo branch that was checked out before we started.
- <b>`made_changes`</b> (bool): Whether making the new template pipeline introduced any changes
- <b>`make_pr`</b> (bool): Whether to try to automatically make a PR on GitHub.com
- <b>`required_config_vars`</b> (list): List of nextflow variables required to make template pipeline
- <b>`gh_username`</b> (str): GitHub username
- <b>`gh_repo`</b> (str): GitHub repository name
- <b>`template_yaml`</b> (str): Path to template.yml file for pipeline creation settings.

<a href="../../../../../../tools/nf_core/sync.py#L58"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline_dir,
    from_branch=None,
    make_pr=False,
    gh_repo=None,
    gh_username=None,
    template_yaml_path=None
)
```

Initialise syncing object

---

<a href="../../../../../../tools/nf_core/sync.py#L194"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `checkout_template_branch`

```python
checkout_template_branch()
```

Try to check out the origin/TEMPLATE in a new TEMPLATE branch. If this fails, try to check out an existing local TEMPLATE branch.

---

<a href="../../../../../../tools/nf_core/sync.py#L407"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `close_open_pr`

```python
close_open_pr(pr)
```

Given a PR API response, add a comment and close.

---

<a href="../../../../../../tools/nf_core/sync.py#L369"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `close_open_template_merge_prs`

```python
close_open_template_merge_prs()
```

Get all template merging branches (starting with 'nf-core-template-merge-') and check for any open PRs from these branches to the self.from_branch If open PRs are found, add a comment and close them

---

<a href="../../../../../../tools/nf_core/sync.py#L262"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `commit_template_changes`

```python
commit_template_changes()
```

If we have any changes with the new template files, make a git commit

---

<a href="../../../../../../tools/nf_core/sync.py#L289"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create_merge_base_branch`

```python
create_merge_base_branch()
```

Create a new branch from the updated TEMPLATE branch This branch will then be used to create the PR

---

<a href="../../../../../../tools/nf_core/sync.py#L209"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `delete_template_branch_files`

```python
delete_template_branch_files()
```

Delete all files in the TEMPLATE branch

---

<a href="../../../../../../tools/nf_core/sync.py#L166"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_wf_config`

```python
get_wf_config()
```

Check out the target branch if requested and fetch the nextflow config. Check that we have the required config variables.

---

<a href="../../../../../../tools/nf_core/sync.py#L146"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `inspect_sync_dir`

```python
inspect_sync_dir()
```

Takes a look at the target directory for syncing. Checks that it's a git repo and makes sure that there are no uncommitted changes.

---

<a href="../../../../../../tools/nf_core/sync.py#L325"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_pull_request`

```python
make_pull_request()
```

Create a pull request to a base branch (default: dev), from a head branch (default: TEMPLATE)

Returns: An instance of class requests.Response

---

<a href="../../../../../../tools/nf_core/sync.py#L228"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_template_pipeline`

```python
make_template_pipeline()
```

Delete all files and make a fresh template using the workflow variables

---

<a href="../../../../../../tools/nf_core/sync.py#L316"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `push_merge_branch`

```python
push_merge_branch()
```

Push the newly created merge branch to the remote repository

---

<a href="../../../../../../tools/nf_core/sync.py#L278"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `push_template_branch`

```python
push_template_branch()
```

If we made any changes, push the TEMPLATE branch to the default remote and try to make a PR. If we don't have the auth token, try to figure out a URL for the PR and print this to the console.

---

<a href="../../../../../../tools/nf_core/sync.py#L442"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `reset_target_dir`

```python
reset_target_dir()
```

Reset the target pipeline directory. Check out the original branch.

---

<a href="../../../../../../tools/nf_core/sync.py#L98"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `sync`

```python
sync()
```

Find workflow attributes, create a new template pipeline on TEMPLATE

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
