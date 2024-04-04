# nf\_core.download

Downloads a nf-core pipeline to the local file system.

### *`exception{:python}`*`nf_core.download.ContainerError(container, registry, address, absolute_URI, out_path, singularity_command, error_msg){:python}`

Bases: `Exception`

A class of errors related to pulling containers with Singularity/Apptainer

#### *`exception{:python}`*`ImageExistsError(error_log){:python}`

Bases: `FileExistsError`

Image already exists in cache/output directory.

#### *`exception{:python}`*`ImageNotFoundError(error_log){:python}`

Bases: `FileNotFoundError`

The image can not be found in the registry

#### *`exception{:python}`*`InvalidTagError(error_log){:python}`

Bases: `AttributeError`

Image and registry are valid, but the (version) tag is not

#### *`exception{:python}`*`OtherError(error_log){:python}`

Bases: `RuntimeError`

Undefined error with the container

#### *`exception{:python}`*`RegistryNotFoundError(error_log){:python}`

Bases: `ConnectionRefusedError`

The specified registry does not resolve to a valid IP address

### *`exception{:python}`*`nf_core.download.DownloadError{:python}`

Bases: `RuntimeError`

A custom exception that is raised when nf-core download encounters a problem that we already took into consideration.
In this case, we do not want to print the traceback, but give the user some concise, helpful feedback instead.

### *`class{:python}`*`nf_core.download.DownloadProgress(*columns: str | ProgressColumn, console: Console | None = None, auto_refresh: bool = True, refresh_per_second: float = 10, speed_estimate_period: float = 30.0, transient: bool = False, redirect_stdout: bool = True, redirect_stderr: bool = True, get_time: Callable[[], float] | None = None, disable: bool = False, expand: bool = False){:python}`

Bases: `Progress`

Custom Progress bar class, allowing us to have two progress
bars with different columns / layouts.

#### `get_renderables(){:python}`

Get a number of renderables for the progress display.

### *`class{:python}`*`nf_core.download.DownloadWorkflow(pipeline=None, revision=None, outdir=None, compress_type=None, force=False, tower=False, download_configuration=None, container_system=None, container_library=None, container_cache_utilisation=None, container_cache_index=None, parallel_downloads=4){:python}`

Bases: `object`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

* **Parameters:**
  * **pipeline** (*str*) – A nf-core pipeline name.
  * **revision** (*List* \*\[\**str* *]*) – The workflow revision to download, like 1.0. Defaults to None.
  * **container** (*bool*) – Flag, if the Singularity container should be downloaded as well. Defaults to False.
  * **tower** (*bool*) – Flag, to customize the download for Nextflow Tower (convert to git bare repo). Defaults to False.
  * **outdir** (*str*) – Path to the local download directory. Defaults to None.

#### `compress_download(){:python}`

Take the downloaded files and make a compressed .tar.gz archive.

#### `download_configs(){:python}`

Downloads the centralised config profiles from nf-core/configs to `self.outdir`.

#### `download_wf_files(revision, wf_sha, download_url){:python}`

Downloads workflow files from GitHub to the `self.outdir`.

#### `download_workflow(){:python}`

Starts a nf-core workflow download.

#### `download_workflow_static(){:python}`

Downloads a nf-core workflow from GitHub to the local file system in a self-contained manner.

#### `download_workflow_tower(location=None){:python}`

Create a bare-cloned git repository of the workflow, so it can be launched with tw launch as <file:/> pipeline

#### `find_container_images(workflow_directory){:python}`

Find container image names for workflow.

Starts by using nextflow config to pull out any process.container
declarations. This works for DSL1. It should return a simple string with resolved logic,
but not always, e.g. not for differentialabundance 1.2.0

Second, we look for DSL2 containers. These can’t be found with
nextflow config at the time of writing, so we scrape the pipeline files.
This returns raw matches that will likely need to be cleaned.

#### `gather_registries(workflow_directory: str){:python}`

Fetch the registries from the pipeline config and CLI arguments and store them in a set.
This is needed to symlink downloaded container images so Nextflow will find them.

#### `get_revision_hash(){:python}`

Find specified revision / branch hash

#### `get_singularity_images(current_revision: str = ''){:python}`

Loop through container names and download Singularity images

#### `prioritize_direct_download(container_list){:python}`

Helper function that takes a list of container images (URLs and Docker URIs),
eliminates all Docker URIs for which also a URL is contained and returns the
cleaned and also deduplicated list.

Conceptually, this works like so:

Everything after the last Slash should be identical, e.g. “scanpy:1.7.2–pyhdfd78af\_0” in
\[‘https://depot.galaxyproject.org/singularity/scanpy:1.7.2–pyhdfd78af\_0’, ‘biocontainers/scanpy:1.7.2–pyhdfd78af\_0’]

re.sub(‘.\*/(.\*)’,’1’,c) will drop everything up to the last slash from c (container\_id)

d.get(k:=re.sub(‘.\*/(.\*)’,’1’,c),’’) assigns the truncated string to k (key) and gets the
corresponding value from the dict if present or else defaults to “”.

If the regex pattern matches, the original container\_id will be assigned to the dict with the k key.
r”^$|(?!^http)” matches an empty string (we didn’t have it in the dict yet and want to keep it in either case) or
any string that does not start with http. Because if our current dict value already starts with http,
we want to keep it and not replace with with whatever we have now (which might be the Docker URI).

A regex that matches http, r”^$|^http” could thus be used to prioritize the Docker URIs over http Downloads

#### `prompt_compression_type(){:python}`

Ask user if we should compress the downloaded files

#### `prompt_config_inclusion(){:python}`

Prompt for inclusion of institutional configurations

#### `prompt_container_download(){:python}`

Prompt whether to download container images or not

#### `prompt_pipeline_name(){:python}`

Prompt for the pipeline name if not set with a flag

#### `prompt_revision(){:python}`

Prompt for pipeline revision / branch
Prompt user for revision tag if ‘–revision’ was not set
If –tower is specified, allow to select multiple revisions
Also the static download allows for multiple revisions, but
we do not prompt this option interactively.

#### `prompt_singularity_cachedir_creation(){:python}`

Prompt about using $NXF\_SINGULARITY\_CACHEDIR if not already set

#### `prompt_singularity_cachedir_remote(){:python}`

Prompt about the index of a remote $NXF\_SINGULARITY\_CACHEDIR

#### `prompt_singularity_cachedir_utilization(){:python}`

Ask if we should *only* use $NXF\_SINGULARITY\_CACHEDIR without copying into target

#### `read_remote_containers(){:python}`

Reads the file specified as index for the remote Singularity cache dir

#### `rectify_raw_container_matches(raw_findings){:python}`

Helper function to rectify the raw extracted container matches into fully qualified container names.
If multiple containers are found, any prefixed with http for direct download is prioritized

Example syntax:

Early DSL2:

```groovy
if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
} else {
    container "quay.io/biocontainers/fastqc:0.11.9--0"
}
```

Later DSL2:

```groovy
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    'biocontainers/fastqc:0.11.9--0' }"
```

Later DSL2, variable is being used:

```groovy
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    "https://depot.galaxyproject.org/singularity/${container_id}" :
    "quay.io/biocontainers/${container_id}" }"

container_id = 'mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0'
```

DSL1 / Special case DSL2:

```groovy
container "nfcore/cellranger:6.0.2"
```

#### `singularity_copy_cache_image(container: str, out_path: str, cache_path: str | None){:python}`

Copy Singularity image from NXF\_SINGULARITY\_CACHEDIR to target folder.

#### `singularity_download_image(container: str, out_path: str, cache_path: str | None, progress:{:python}`[DownloadProgress](#nf_core.download.DownloadProgress))

Download a singularity image from the web.

Use native Python to download the file.

* **Parameters:**
  * **container** (*str*) – A pipeline’s container name. Usually it is of similar format
    to `https://depot.galaxyproject.org/singularity/name:version`
  * **out\_path** (*str*) – The final target output path
  * **cache\_path** (*str* *,* *None*) – The NXF\_SINGULARITY\_CACHEDIR path if set, None if not
  * **progress** (*Progress*) – Rich progress bar instance to add tasks to.

#### `singularity_image_filenames(container: str){:python}`

Check Singularity cache for image, copy to destination folder if found.

* **Parameters:**
  **container** (*str*) – A pipeline’s container name. Can be direct download URL
  or a Docker Hub repository ID.
* **Returns:**
  Returns a tuple of (out\_path, cache\_path).
  : out\_path is the final target output path. it may point to the NXF\_SINGULARITY\_CACHEDIR, if cache utilisation was set to ‘amend’.
  If cache utilisation was set to ‘copy’, it will point to the target folder, a subdirectory of the output directory. In the latter case,
  cache\_path may either be None (image is not yet cached locally) or point to the image in the NXF\_SINGULARITY\_CACHEDIR, so it will not be
  downloaded from the web again, but directly copied from there. See get\_singularity\_images() for implementation.
* **Return type:**
  tuple (str, str)

#### `singularity_pull_image(container: str, out_path: str, cache_path: str | None, library: List[str], progress:{:python}`[DownloadProgress](#nf_core.download.DownloadProgress))

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

* **Parameters:**
  * **container** (*str*) – A pipeline’s container name. Usually it is of similar format
    to `nfcore/name:version`.
  * **library** (*list* *of* *str*) – A list of libraries to try for pulling the image.
* **Raises:**
  **Various exceptions possible from subprocess execution** **of** **Singularity.** –

#### `symlink_singularity_images(image_out_path: str){:python}`

Create a symlink for each registry in the registry set that points to the image.
We have dropped the explicit registries from the modules in favor of the configurable registries.
Unfortunately, Nextflow still expects the registry to be part of the file name, so a symlink is needed.

The base image, e.g. ./nf-core-gatk-4.4.0.0.img will thus be symlinked as for example ./quay.io-nf-core-gatk-4.4.0.0.img
by prepending all registries in self.registry\_set to the image name.

Unfortunately, out output image name may contain a registry definition (Singularity image pulled from depot.galaxyproject.org
or older pipeline version, where the docker registry was part of the image name in the modules). Hence, it must be stripped
before to ensure that it is really the base name.

#### `wf_use_local_configs(revision_dirname){:python}`

Edit the downloaded nextflow.config file to use the local config files

### *`class{:python}`*`nf_core.download.WorkflowRepo(remote_url, revision, commit, location=None, hide_progress=False, in_cache=True){:python}`

Bases: `SyncedRepo`

An object to store details about a locally cached workflow repository.

Important Attributes:
: fullname: The full name of the repository, `nf-core/{self.pipelinename}`.
local\_repo\_dir (str): The local directory, where the workflow is cloned into. Defaults to `$HOME/.cache/nf-core/nf-core/{self.pipeline}`.

#### `access(){:python}`

#### `bare_clone(destination){:python}`

#### `checkout(commit){:python}`

Checks out the repository at the requested commit

* **Parameters:**
  **commit** (*str*) – Git SHA of the commit

#### `get_remote_branches(remote_url){:python}`

Get all branches from a remote repository

* **Parameters:**
  **remote\_url** (*str*) – The git url to the remote repository
* **Returns:**
  All branches found in the remote
* **Return type:**
  (set\[str])

#### `retry_setup_local_repo(skip_confirm=False){:python}`

#### `setup_local_repo(remote, location=None, in_cache=True){:python}`

Sets up the local git repository. If the repository has been cloned previously, it
returns a git.Repo object of that clone. Otherwise it tries to clone the repository from
the provided remote URL and returns a git.Repo of the new clone.

* **Parameters:**
  * **remote** (*str*) – git url of remote
  * **location** (*Path*) – location where the clone should be created/cached.
  * **in\_cache** (*bool* *,* *optional*) – Whether to clone the repository from the cache. Defaults to False.

Sets self.repo

#### `tidy_tags_and_branches(){:python}`

Function to delete all tags and branches that are not of interest to the downloader.
This allows a clutter-free experience in Tower. The untagged commits are evidently still available.

However, due to local caching, the downloader might also want access to revisions that had been deleted before.
In that case, don’t bother with re-adding the tags and rather download  anew from Github.
