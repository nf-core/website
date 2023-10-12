<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/sync.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.sync`

Synchronise a pipeline TEMPLATE branch with the template.

---

<a href="../../../../../../tools/nf_core/sync.py#L370"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `sync_all_pipelines`

```python
sync_all_pipelines(gh_username=None, gh_auth_token=None)
```

Sync all nf-core pipelines

---

<a href="../../../../../../tools/nf_core/sync.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SyncException`

Exception raised when there was an error with TEMPLATE branch synchronisation

---

<a href="../../../../../../tools/nf_core/sync.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PullRequestException`

Exception raised when there was an error creating a Pull-Request on GitHub.com

---

<a href="../../../../../../tools/nf_core/sync.py#L27"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineSync`

Object to hold syncing information and results.

**Args:**

- <b>`pipeline_dir`</b> (str): The path to the Nextflow pipeline root directory
- <b>`make_template_branch`</b> (bool): Set this to `True` to create a `TEMPLATE` branch if it is not found
- <b>`from_branch`</b> (str): The branch to use to fetch config vars. If not set, will use current active branch
- <b>`make_pr`</b> (bool): Set this to `True` to create a GitHub pull-request with the changes
- <b>`gh_username`</b> (str): GitHub username
- <b>`gh_repo`</b> (str): GitHub repository name
- <b>`gh_auth_token`</b> (str): Authorisation token used to make PR with GitHub API

**Attributes:**

- <b>`pipeline_dir`</b> (str): Path to target pipeline directory
- <b>`from_branch`</b> (str): Repo branch to use when collecting workflow variables. Default: active branch.
- <b>`make_template_branch`</b> (bool): Whether to try to create TEMPLATE branch if not found
- <b>`orphan_branch`</b> (bool): Whether an orphan branch was made when creating TEMPLATE
- <b>`made_changes`</b> (bool): Whether making the new template pipeline introduced any changes
- <b>`make_pr`</b> (bool): Whether to try to automatically make a PR on GitHub.com
- <b>`required_config_vars`</b> (list): List of nextflow variables required to make template pipeline
- <b>`gh_username`</b> (str): GitHub username
- <b>`gh_repo`</b> (str): GitHub repository name
- <b>`gh_auth_token`</b> (str): Authorisation token used to make PR with GitHub API

<a href="../../../../../../tools/nf_core/sync.py#L52"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline_dir,
    make_template_branch=False,
    from_branch=None,
    make_pr=False,
    gh_username=None,
    gh_repo=None,
    gh_auth_token=None
)
```

Initialise syncing object

---

<a href="../../../../../../tools/nf_core/sync.py#L176"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `checkout_template_branch`

```python
checkout_template_branch()
```

Try to check out the TEMPLATE branch. If it fails, try origin/TEMPLATE. If it still fails and --make-template-branch was given, create it as an orphan branch.

---

<a href="../../../../../../tools/nf_core/sync.py#L250"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `commit_template_changes`

```python
commit_template_changes()
```

If we have any changes with the new template files, make a git commit

---

<a href="../../../../../../tools/nf_core/sync.py#L134"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_wf_config`

```python
get_wf_config()
```

Check out the target branch if requested and fetch the nextflow config. Check that we have the required config variables.

---

<a href="../../../../../../tools/nf_core/sync.py#L350"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `git_merge_help`

```python
git_merge_help()
```

Print a command line help message with instructions on how to merge changes

---

<a href="../../../../../../tools/nf_core/sync.py#L116"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `inspect_sync_dir`

```python
inspect_sync_dir()
```

Takes a look at the target directory for syncing. Checks that it's a git repo and makes sure that there are no uncommitted changes.

---

<a href="../../../../../../tools/nf_core/sync.py#L285"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_pull_request`

```python
make_pull_request()
```

Create a pull request to a base branch (default: dev), from a head branch (default: TEMPLATE)

Returns: An instance of class requests.Response

---

<a href="../../../../../../tools/nf_core/sync.py#L210"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_template_pipeline`

```python
make_template_pipeline()
```

Delete all files and make a fresh template using the workflow variables

---

<a href="../../../../../../tools/nf_core/sync.py#L265"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `push_template_branch`

```python
push_template_branch()
```

If we made any changes, push the TEMPLATE branch to the default remote and try to make a PR. If we don't have the auth token, try to figure out a URL for the PR and print this to the console.

---

<a href="../../../../../../tools/nf_core/sync.py#L339"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `reset_target_dir`

```python
reset_target_dir()
```

Reset the target pipeline directory. Check out the original branch.

---

<a href="../../../../../../tools/nf_core/sync.py#L74"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `sync`

```python
sync()
```

Find workflow attributes, create a new template pipeline on TEMPLATE

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
