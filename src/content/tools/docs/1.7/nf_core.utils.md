<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.utils`

Common utility functions for the nf-core python package.

---

<a href="../../../../../../tools/nf_core/utils.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetch_wf_config`

```python
fetch_wf_config(wf_path, wf=None)
```

Uses Nextflow to retrieve the the configuration variables from a Nextflow workflow.

**Args:**

- <b>`wf_path`</b> (str): Nextflow workflow file system path.

**Returns:**

- <b>`dict`</b>: Workflow configuration settings.

---

<a href="../../../../../../tools/nf_core/utils.py#L71"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_requests_cachedir`

```python
setup_requests_cachedir()
```

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user's home directory under a .nfcore_cache subdir.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
