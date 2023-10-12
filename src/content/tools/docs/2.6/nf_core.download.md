<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.download`

Downloads a nf-core pipeline to the local file system.

---

<a href="../../../../../../tools/nf_core/download.py#L34"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_renderables`

```python
get_renderables()
```

---

<a href="../../../../../../tools/nf_core/download.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadWorkflow`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

**Args:**

- <b>`pipeline`</b> (str): A nf-core pipeline name.
- <b>`revision`</b> (str): The workflow revision to download, like `1.0`. Defaults to None.
- <b>`singularity`</b> (bool): Flag, if the Singularity container should be downloaded as well. Defaults to False.
- <b>`outdir`</b> (str): Path to the local download directory. Defaults to None.

<a href="../../../../../../tools/nf_core/download.py#L80"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline=None,
    revision=None,
    outdir=None,
    compress_type=None,
    force=False,
    container=None,
    singularity_cache_only=False,
    parallel_downloads=4
)
```

---

<a href="../../../../../../tools/nf_core/download.py#L771"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compress_download`

```python
compress_download()
```

Take the downloaded files and make a compressed .tar.gz archive.

---

<a href="../../../../../../tools/nf_core/download.py#L365"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_configs`

```python
download_configs()
```

Downloads the centralised config profiles from nf-core/configs to :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L347"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_wf_files`

```python
download_wf_files()
```

Downloads workflow files from GitHub to the :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow`

```python
download_workflow()
```

Starts a nf-core workflow download.

---

<a href="../../../../../../tools/nf_core/download.py#L412"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_container_images`

```python
find_container_images()
```

Find container image names for workflow.

Starts by using `nextflow config` to pull out any process.container declarations. This works for DSL1. It should return a simple string with resolved logic.

Second, we look for DSL2 containers. These can't be found with `nextflow config` at the time of writing, so we scrape the pipeline files. This returns raw source code that will likely need to be cleaned.

If multiple containers are found, prioritise any prefixed with http for direct download.

Example syntax:

Early DSL2: if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) { container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" } else { container "quay.io/biocontainers/fastqc:0.11.9--0" }

Later DSL2: container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' : 'quay.io/biocontainers/fastqc:0.11.9--0' }"

DSL1 / Special case DSL2: container "nfcore/cellranger:6.0.2"

---

<a href="../../../../../../tools/nf_core/download.py#L202"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_revision_hash`

```python
get_revision_hash()
```

Find specified revision / branch hash

---

<a href="../../../../../../tools/nf_core/download.py#L495"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_singularity_images`

```python
get_singularity_images()
```

Loop through container names and download Singularity images

---

<a href="../../../../../../tools/nf_core/download.py#L322"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_compression_type`

```python
prompt_compression_type()
```

Ask user if we should compress the downloaded files

---

<a href="../../../../../../tools/nf_core/download.py#L233"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_container_download`

```python
prompt_container_download()
```

Prompt whether to download container images or not

---

<a href="../../../../../../tools/nf_core/download.py#L189"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_pipeline_name`

```python
prompt_pipeline_name()
```

Prompt for the pipeline name if not set with a flag

---

<a href="../../../../../../tools/nf_core/download.py#L196"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_revision`

```python
prompt_revision()
```

Prompt for pipeline revision / branch

---

<a href="../../../../../../tools/nf_core/download.py#L302"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_singularity_cachedir_only`

```python
prompt_singularity_cachedir_only()
```

Ask if we should _only_ use $NXF_SINGULARITY_CACHEDIR without copying into target

---

<a href="../../../../../../tools/nf_core/download.py#L244"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_use_singularity_cachedir`

```python
prompt_use_singularity_cachedir()
```

Prompt about using $NXF_SINGULARITY_CACHEDIR if not already set

---

<a href="../../../../../../tools/nf_core/download.py#L643"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_copy_cache_image`

```python
singularity_copy_cache_image(container, out_path, cache_path)
```

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

---

<a href="../../../../../../tools/nf_core/download.py#L650"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L598"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L720"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L384"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `wf_use_local_configs`

```python
wf_use_local_configs()
```

Edit the downloaded nextflow.config file to use the local config files

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
