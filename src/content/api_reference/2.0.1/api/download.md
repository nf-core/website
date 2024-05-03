# nf_core.download

Downloads a nf-core pipeline to the local file system.

### _`class{:python}`_`nf_core.download.DownloadProgress(*columns: str | ProgressColumn, console: Console | None = None, auto_refresh: bool = True, refresh_per_second: float = 10, speed_estimate_period: float = 30.0, transient: bool = False, redirect_stdout: bool = True, redirect_stderr: bool = True, get_time: Callable[[], float] | None = None, disable: bool = False, expand: bool = False){:python}`

Bases: `Progress`

Custom Progress bar class, allowing us to have two progress
bars with different columns / layouts.

#### `get_renderables(){:python}`

Get a number of renderables for the progress display.

### _`class{:python}`_`nf_core.download.DownloadWorkflow(pipeline=None, release=None, outdir=None, compress_type=None, force=False, container=None, singularity_cache_only=False, parallel_downloads=4){:python}`

Bases: `object`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

- **Parameters:**
  - **pipeline** (_str_) – A nf-core pipeline name.
  - **release** (_str_) – The workflow release version to download, like 1.0. Defaults to None.
  - **singularity** (_bool_) – Flag, if the Singularity container should be downloaded as well. Defaults to False.
  - **outdir** (_str_) – Path to the local download directory. Defaults to None.

#### `compress_download(){:python}`

Take the downloaded files and make a compressed .tar.gz archive.

#### `download_configs(){:python}`

Downloads the centralised config profiles from nf-core/configs to `self.outdir`.

#### `download_wf_files(){:python}`

Downloads workflow files from GitHub to the `self.outdir`.

#### `download_workflow(){:python}`

Starts a nf-core workflow download.

#### `find_container_images(){:python}`

Find container image names for workflow.

Starts by using nextflow config to pull out any process.container
declarations. This works for DSL1.

Second, we look for DSL2 containers. These can’t be found with
nextflow config at the time of writing, so we scrape the pipeline files.

#### `get_release_hash(){:python}`

Find specified release / branch hash

#### `get_singularity_images(){:python}`

Loop through container names and download Singularity images

#### `prompt_compression_type(){:python}`

Ask user if we should compress the downloaded files

#### `prompt_container_download(){:python}`

Prompt whether to download container images or not

#### `prompt_pipeline_name(){:python}`

Prompt for the pipeline name if not set with a flag

#### `prompt_release(){:python}`

Prompt for pipeline release / branch

#### `prompt_singularity_cachedir_only(){:python}`

Ask if we should _only_ use $NXF_SINGULARITY_CACHEDIR without copying into target

#### `prompt_use_singularity_cachedir(){:python}`

Prompt about using $NXF_SINGULARITY_CACHEDIR if not already set

#### `singularity_copy_cache_image(container, out_path, cache_path){:python}`

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

#### `singularity_download_image(container, out_path, cache_path, progress){:python}`

Download a singularity image from the web.

Use native Python to download the file.

- **Parameters:**
  - **container** (_str_) – A pipeline’s container name. Usually it is of similar format
    to `https://depot.galaxyproject.org/singularity/name:version`
  - **out_path** (_str_) – The final target output path
  - **cache_path** (_str_\*,\* _None_) – The NXF_SINGULARITY_CACHEDIR path if set, None if not
  - **progress** (_Progress_) – Rich progress bar instance to add tasks to.

#### `singularity_image_filenames(container){:python}`

Check Singularity cache for image, copy to destination folder if found.

- **Parameters:**
  **container** (_str_) – A pipeline’s container name. Can be direct download URL
  or a Docker Hub repository ID.
- **Returns:**
  Returns True if we have the image in the target location.
  : Returns a download path if not.
- **Return type:**
  results (bool, str)

#### `singularity_pull_image(container, out_path, cache_path, progress){:python}`

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

- **Parameters:**
  **container** (_str_) – A pipeline’s container name. Usually it is of similar format
  to `nfcore/name:version`.
- **Raises:**
  **Various exceptions possible from subprocess execution** **of** **Singularity.** –

#### `validate_md5(fname, expected=None){:python}`

Calculates the md5sum for a file on the disk and validate with expected.

- **Parameters:**
  - **fname** (_str_) – Path to a local file.
  - **expected** (_str_) – The expected md5sum.
- **Raises:**
  **IOError**\*\*,\*\* **if the md5sum does not match the remote sum.** –

#### `wf_use_local_configs(){:python}`

Edit the downloaded nextflow.config file to use the local config files
