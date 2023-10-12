<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.download`

Downloads a nf-core pipeline to the local file system.

---

<a href="../../../../../../tools/nf_core/download.py#L29"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadProgress`

Custom Progress bar class, allowing us to have two progress bars with different columns / layouts.

---

#### <kbd>property</kbd> console

---

#### <kbd>property</kbd> finished

Check if all tasks have been completed.

---

#### <kbd>property</kbd> task_ids

A list of task IDs.

---

#### <kbd>property</kbd> tasks

Get a list of Task instances.

---

<a href="../../../../../../tools/nf_core/download.py#L34"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_renderables`

```python
get_renderables()
```

---

<a href="../../../../../../tools/nf_core/download.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadWorkflow`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

**Args:**

- <b>`pipeline`</b> (str): A nf-core pipeline name.
- <b>`release`</b> (str): The workflow release version to download, like `1.0`. Defaults to None.
- <b>`singularity`</b> (bool): Flag, if the Singularity container should be downloaded as well. Defaults to False.
- <b>`outdir`</b> (str): Path to the local download directory. Defaults to None.

<a href="../../../../../../tools/nf_core/download.py#L75"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline,
    release=None,
    outdir=None,
    compress_type='tar.gz',
    force=False,
    singularity=False,
    singularity_cache_only=False,
    parallel_downloads=4
)
```

---

<a href="../../../../../../tools/nf_core/download.py#L643"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compress_download`

```python
compress_download()
```

Take the downloaded files and make a compressed .tar.gz archive.

---

<a href="../../../../../../tools/nf_core/download.py#L266"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_configs`

```python
download_configs()
```

Downloads the centralised config profiles from nf-core/configs to :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L248"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_wf_files`

```python
download_wf_files()
```

Downloads workflow files from GitHub to the :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L109"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow`

```python
download_workflow()
```

Starts a nf-core workflow download.

---

<a href="../../../../../../tools/nf_core/download.py#L169"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L311"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_container_images`

```python
find_container_images()
```

Find container image names for workflow.

Starts by using `nextflow config` to pull out any process.container declarations. This works for DSL1.

Second, we look for DSL2 containers. These can't be found with `nextflow config` at the time of writing, so we scrape the pipeline files.

---

<a href="../../../../../../tools/nf_core/download.py#L358"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_singularity_images`

```python
get_singularity_images()
```

Loop through container names and download Singularity images

---

<a href="../../../../../../tools/nf_core/download.py#L513"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_copy_cache_image`

```python
singularity_copy_cache_image(container, out_path, cache_path)
```

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

---

<a href="../../../../../../tools/nf_core/download.py#L520"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_download_image`

```python
singularity_download_image(container, out_path, cache_path, progress)
```

Download a singularity image from the web.

Use native Python to download the file.

**Args:**

- <b>`container`</b> (str): A pipeline's container name. Usually it is of similar format
- <b>` to ``https `</b>: //depot.galaxyproject.org/singularity/name:version``
- <b>`out_path`</b> (str): The final target output path
- <b>`cache_path`</b> (str, None): The NXF_SINGULARITY_CACHEDIR path if set, None if not
- <b>`progress`</b> (Progress): Rich progress bar instance to add tasks to.

---

<a href="../../../../../../tools/nf_core/download.py#L466"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_image_filenames`

```python
singularity_image_filenames(container)
```

Check Singularity cache for image, copy to destination folder if found.

**Args:**

- <b>`container`</b> (str): A pipeline's container name. Can be direct download URL or a Docker Hub repository ID.

**Returns:**

- <b>`results`</b> (bool, str): Returns True if we have the image in the target location. Returns a download path if not.

---

<a href="../../../../../../tools/nf_core/download.py#L590"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_pull_image`

```python
singularity_pull_image(container, out_path, cache_path, progress)
```

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

**Args:**

- <b>`container`</b> (str): A pipeline's container name. Usually it is of similar format
- <b>` to ``nfcore/name `</b>: version``.

**Raises:**
Various exceptions possible from `subprocess` execution of Singularity.

---

<a href="../../../../../../tools/nf_core/download.py#L674"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L285"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `wf_use_local_configs`

```python
wf_use_local_configs()
```

Edit the downloaded nextflow.config file to use the local config files

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
