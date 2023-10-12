<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.download`

Downloads a nf-core pipeline to the local file system.

## **Global Variables**

- **NFCORE_CACHE_DIR**
- **NFCORE_DIR**

---

<a href="../../../../../../tools/nf_core/download.py#L42"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadError`

A custom exception that is raised when nf-core download encounters a problem that we already took into consideration. In this case, we do not want to print the traceback, but give the user some concise, helpful feedback instead.

---

<a href="../../../../../../tools/nf_core/download.py#L47"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L52"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_renderables`

```python
get_renderables()
```

---

<a href="../../../../../../tools/nf_core/download.py#L81"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadWorkflow`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

**Args:**

- <b>`pipeline`</b> (str): A nf-core pipeline name.
- <b>`revision`</b> (List[str]): The workflow revision to download, like `1.0`. Defaults to None.
- <b>`container`</b> (bool): Flag, if the Singularity container should be downloaded as well. Defaults to False.
- <b>`tower`</b> (bool): Flag, to customize the download for Nextflow Tower (convert to git bare repo). Defaults to False.
- <b>`outdir`</b> (str): Path to the local download directory. Defaults to None.

<a href="../../../../../../tools/nf_core/download.py#L94"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    pipeline=None,
    revision=None,
    outdir=None,
    compress_type=None,
    force=False,
    tower=False,
    download_configuration=None,
    container_system=None,
    container_library=None,
    container_cache_utilisation=None,
    container_cache_index=None,
    parallel_downloads=4
)
```

---

<a href="../../../../../../tools/nf_core/download.py#L1252"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compress_download`

```python
compress_download()
```

Take the downloaded files and make a compressed .tar.gz archive.

---

<a href="../../../../../../tools/nf_core/download.py#L616"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_configs`

```python
download_configs()
```

Downloads the centralised config profiles from nf-core/configs to :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L590"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_wf_files`

```python
download_wf_files(revision, wf_sha, download_url)
```

Downloads workflow files from GitHub to the :attr:`self.outdir`.

---

<a href="../../../../../../tools/nf_core/download.py#L152"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow`

```python
download_workflow()
```

Starts a nf-core workflow download.

---

<a href="../../../../../../tools/nf_core/download.py#L236"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow_static`

```python
download_workflow_static()
```

Downloads a nf-core workflow from GitHub to the local file system in a self-contained manner.

---

<a href="../../../../../../tools/nf_core/download.py#L270"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow_tower`

```python
download_workflow_tower(location=None)
```

Create a bare-cloned git repository of the workflow, so it can be launched with `tw launch` as file:/ pipeline

---

<a href="../../../../../../tools/nf_core/download.py#L663"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_container_images`

```python
find_container_images(workflow_directory)
```

Find container image names for workflow.

Starts by using `nextflow config` to pull out any process.container declarations. This works for DSL1. It should return a simple string with resolved logic, but not always, e.g. not for differentialabundance 1.2.0

Second, we look for DSL2 containers. These can't be found with `nextflow config` at the time of writing, so we scrape the pipeline files. This returns raw matches that will likely need to be cleaned.

---

<a href="../../../../../../tools/nf_core/download.py#L347"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_revision_hash`

```python
get_revision_hash()
```

Find specified revision / branch hash

---

<a href="../../../../../../tools/nf_core/download.py#L917"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_singularity_images`

```python
get_singularity_images(current_revision='')
```

Loop through container names and download Singularity images

---

<a href="../../../../../../tools/nf_core/download.py#L565"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_compression_type`

```python
prompt_compression_type()
```

Ask user if we should compress the downloaded files

---

<a href="../../../../../../tools/nf_core/download.py#L387"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_config_inclusion`

```python
prompt_config_inclusion()
```

Prompt for inclusion of institutional configurations

---

<a href="../../../../../../tools/nf_core/download.py#L398"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_container_download`

```python
prompt_container_download()
```

Prompt whether to download container images or not

---

<a href="../../../../../../tools/nf_core/download.py#L306"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_pipeline_name`

```python
prompt_pipeline_name()
```

Prompt for the pipeline name if not set with a flag

---

<a href="../../../../../../tools/nf_core/download.py#L313"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_revision`

```python
prompt_revision()
```

Prompt for pipeline revision / branch Prompt user for revision tag if '--revision' was not set If --tower is specified, allow to select multiple revisions Also the static download allows for multiple revisions, but we do not prompt this option interactively.

---

<a href="../../../../../../tools/nf_core/download.py#L409"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_singularity_cachedir_creation`

```python
prompt_singularity_cachedir_creation()
```

Prompt about using $NXF_SINGULARITY_CACHEDIR if not already set

---

<a href="../../../../../../tools/nf_core/download.py#L505"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_singularity_cachedir_remote`

```python
prompt_singularity_cachedir_remote()
```

Prompt about the index of a remote $NXF_SINGULARITY_CACHEDIR

---

<a href="../../../../../../tools/nf_core/download.py#L486"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_singularity_cachedir_utilization`

```python
prompt_singularity_cachedir_utilization()
```

Ask if we should _only_ use $NXF_SINGULARITY_CACHEDIR without copying into target

---

<a href="../../../../../../tools/nf_core/download.py#L534"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `read_remote_containers`

```python
read_remote_containers()
```

Reads the file specified as index for the remote Singularity cache dir

---

<a href="../../../../../../tools/nf_core/download.py#L754"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `rectify_raw_container_matches`

```python
rectify_raw_container_matches(raw_findings)
```

Helper function to rectify the raw extracted container matches into fully qualified container names. If multiple containers are found, any prefixed with http for direct download is prioritized

Example syntax:

Early DSL2: if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) { container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" } else { container "quay.io/biocontainers/fastqc:0.11.9--0" }

Later DSL2: container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' : 'biocontainers/fastqc:0.11.9--0' }"

Later DSL2, variable is being used: container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  "https://depot.galaxyproject.org/singularity/${container_id}" : "quay.io/biocontainers/${container_id}" }"

container_id = 'mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0'

DSL1 / Special case DSL2: container "nfcore/cellranger:6.0.2"

---

<a href="../../../../../../tools/nf_core/download.py#L1113"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_copy_cache_image`

```python
singularity_copy_cache_image(container, out_path, cache_path)
```

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

---

<a href="../../../../../../tools/nf_core/download.py#L1120"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L1068"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/download.py#L1190"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `singularity_pull_image`

```python
singularity_pull_image(container, out_path, cache_path, library, progress)
```

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

**Args:**

- <b>`container`</b> (str): A pipeline's container name. Usually it is of similar format
- <b>` to ``nfcore/name `</b>: version``.
- <b>`library`</b> (list of str): A list of libraries to try for pulling the image.

**Raises:**
Various exceptions possible from `subprocess` execution of Singularity.

---

<a href="../../../../../../tools/nf_core/download.py#L635"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `wf_use_local_configs`

```python
wf_use_local_configs(revision_dirname)
```

Edit the downloaded nextflow.config file to use the local config files

---

<a href="../../../../../../tools/nf_core/download.py#L1284"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `WorkflowRepo`

An object to store details about a locally cached workflow repository.

Important Attributes: fullname: The full name of the repository, `nf-core/{self.pipelinename}`. local_repo_dir (str): The local directory, where the workflow is cloned into. Defaults to `$HOME/.cache/nf-core/nf-core/{self.pipeline}`.

<a href="../../../../../../tools/nf_core/download.py#L1294"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    remote_url,
    revision,
    commit,
    location=None,
    hide_progress=False,
    in_cache=True
)
```

Initializes the object and clones the workflows git repository if it is not already present

**Args:**

- <b>`remote_url`</b> (str): The URL of the remote repository. Defaults to None.
- <b>`self.revision (list of str)`</b>: The revisions to include. A list of strings.
- <b>`commits`</b> (dict of str): The checksums to linked with the revisions.
- <b>`no_pull`</b> (bool, optional): Whether to skip the pull step. Defaults to False.
- <b>`hide_progress`</b> (bool, optional): Whether to hide the progress bar. Defaults to False.
- <b>`in_cache`</b> (bool, optional): Whether to clone the repository from the cache. Defaults to False.

---

<a href="../../../../../../tools/nf_core/download.py#L1340"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `access`

```python
access()
```

---

<a href="../../../../../../tools/nf_core/download.py#L1505"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `bare_clone`

```python
bare_clone(destination)
```

---

<a href="../../../../../../tools/nf_core/download.py#L1346"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `checkout`

```python
checkout(commit)
```

---

<a href="../../../../../../tools/nf_core/download.py#L1349"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_branches`

```python
get_remote_branches(remote_url)
```

---

<a href="../../../../../../tools/nf_core/download.py#L1352"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `retry_setup_local_repo`

```python
retry_setup_local_repo(skip_confirm=False)
```

---

<a href="../../../../../../tools/nf_core/download.py#L1370"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_local_repo`

```python
setup_local_repo(remote, location=None, in_cache=True)
```

Sets up the local git repository. If the repository has been cloned previously, it returns a git.Repo object of that clone. Otherwise it tries to clone the repository from the provided remote URL and returns a git.Repo of the new clone.

**Args:**

- <b>`remote`</b> (str): git url of remote
- <b>`location`</b> (Path): location where the clone should be created/cached.
- <b>`in_cache`</b> (bool, optional): Whether to clone the repository from the cache. Defaults to False. Sets self.repo

---

<a href="../../../../../../tools/nf_core/download.py#L1430"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `tidy_tags_and_branches`

```python
tidy_tags_and_branches()
```

Function to delete all tags and branches that are not of interest to the downloader. This allows a clutter-free experience in Tower. The untagged commits are evidently still available.

However, due to local caching, the downloader might also want access to revisions that had been deleted before. In that case, don't bother with re-adding the tags and rather download anew from Github.

---

<a href="../../../../../../tools/nf_core/download.py#L1521"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ContainerError`

A class of errors related to pulling containers with Singularity/Apptainer

<a href="../../../../../../tools/nf_core/download.py#L1524"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(container, registry, address, out_path, singularity_command, error_msg)
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
