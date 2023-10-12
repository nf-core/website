<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.utils`

Common utility functions for the nf-core python package.

---

<a href="../../../../../../tools/nf_core/utils.py#L23"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_if_outdated`

```python
check_if_outdated(
    current_version=None,
    remote_version=None,
    source_url='https://nf-co.re/tools_version'
)
```

Check if the current version of nf-core is outdated

---

<a href="../../../../../../tools/nf_core/utils.py#L46"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `rich_force_colors`

```python
rich_force_colors()
```

Check if any environment variables are set to force Rich to use coloured output

---

<a href="../../../../../../tools/nf_core/utils.py#L55"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetch_wf_config`

```python
fetch_wf_config(wf_path)
```

Uses Nextflow to retrieve the the configuration variables from a Nextflow workflow.

**Args:**

- <b>`wf_path`</b> (str): Nextflow workflow file system path.

**Returns:**

- <b>`dict`</b>: Workflow configuration settings.

---

<a href="../../../../../../tools/nf_core/utils.py#L140"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_requests_cachedir`

```python
setup_requests_cachedir()
```

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user's home directory under a .nfcore_cache subdir.

---

<a href="../../../../../../tools/nf_core/utils.py#L160"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `wait_cli_function`

```python
wait_cli_function(poll_func, poll_every=20)
```

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

**Arguments:**

- <b>`poll_func`</b> (function): Function to call
- <b>`poll_every`</b> (int): How many tenths of a second to wait between function calls. Default: 20.

**Returns:**
None. Just sits in an infite loop until the function returns True.

---

<a href="../../../../../../tools/nf_core/utils.py#L202"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `poll_nfcore_web_api`

```python
poll_nfcore_web_api(api_url, post_data=None)
```

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level 'status' key.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
