<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.download`

Downloads a nf-core pipeline to the local file system.

---

<a href="../../../../../../tools/nf_core/download.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadWorkflow`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

**Args:**

- <b>`pipeline`</b> (str): A nf-core pipeline name.
- <b>`release`</b> (str): The workflow release version to download, like `1.0`. Defaults to None.
- <b>`singularity`</b> (bool): Flag, if the Singularity container should be downloaded as well. Defaults to False.
- <b>`outdir`</b> (str): Path to the local download directory. Defaults to None.

<a href="../../../../../../tools/nf_core/download.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline,
    release=None,
    singularity=False,
    outdir=None,
    compress_type='tar.gz'
)
```

---

<a href="../../../../../../tools/nf_core/download.py#L286"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compress_download`

```python
compress_download()
```

Take the downloaded files and make a compressed .tar.gz archive.

---

<a href="../../../../../../tools/nf_core/download.py#L206"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_configs`

```python
download_configs()
```

Downloads the centralised config profiles from nf-core/configs to :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L187"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_wf_files`

```python
download_wf_files()
```

Downloads workflow files from GitHub to the :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L49"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow`

```python
download_workflow()
```

Starts a nf-core workflow download.

---

<a href="../../../../../../tools/nf_core/download.py#L114"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `fetch_workflow_details`

```python
fetch_workflow_details(wfs)
```

Fetches details of a nf-core workflow to download.

**Args:**

- <b>`wfs`</b> (nf_core.list.Workflows): A nf_core.list.Workflows object

**Raises:**
LockupError, if the pipeline can not be found.

---

<a href="../../../../../../tools/nf_core/download.py#L246"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_container_images`

```python
find_container_images()
```

Find container image names for workflow

---

<a href="../../../../../../tools/nf_core/download.py#L258"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `pull_singularity_image`

```python
pull_singularity_image(container)
```

Uses a local installation of singularity to pull an image from Docker Hub.

**Args:**

- <b>`container`</b> (str): A pipeline's container name. Usually it is of similar format
- <b>`to `nfcore/name`</b>: dev`.

**Raises:**
Various exceptions possible from `subprocess` execution of Singularity.

---

<a href="../../../../../../tools/nf_core/download.py#L319"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_md5`

```python
validate_md5(fname, expected=None)
```

Calculates the md5sum for a file on the disk and validate with expected.

**Args:**

- <b>`fname`</b> (str): Path to a local file.
- <b>`expected`</b> (str): The expected md5sum.

**Raises:**
IOError, if the md5sum does not match the remote sum.

---

<a href="../../../../../../tools/nf_core/download.py#L226"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `wf_use_local_configs`

```python
wf_use_local_configs()
```

Edit the downloaded nextflow.config file to use the local config files

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
