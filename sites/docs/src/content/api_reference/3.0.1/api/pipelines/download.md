# nf_core.download

Downloads a nf-core pipeline to the local file system.

### _`exception{:python}`_`nf_core.pipelines.download.ContainerError(container, registry, address, absolute_URI, out_path, singularity_command, error_msg){:python}`

Bases: `Exception`

A class of errors related to pulling containers with Singularity/Apptainer

#### _`exception{:python}`_`ImageExistsError(error_log){:python}`

Bases: `FileExistsError`

Image already exists in cache/output directory.

#### _`exception{:python}`_`ImageNotFoundError(error_log){:python}`

Bases: `FileNotFoundError`

The image can not be found in the registry

#### _`exception{:python}`_`InvalidTagError(error_log){:python}`

Bases: `AttributeError`

Image and registry are valid, but the (version) tag is not

#### _`exception{:python}`_`OtherError(error_log){:python}`

Bases: `RuntimeError`

Undefined error with the container

#### _`exception{:python}`_`RegistryNotFoundError(error_log){:python}`

Bases: `ConnectionRefusedError`

The specified registry does not resolve to a valid IP address

### _`exception{:python}`_`nf_core.pipelines.download.DownloadError{:python}`

Bases: `RuntimeError`

A custom exception that is raised when nf-core pipelines download encounters a problem that we already took into consideration.
In this case, we do not want to print the traceback, but give the user some concise, helpful feedback instead.

### _`class{:python}`_`nf_core.pipelines.download.DownloadProgress(*columns: str | ProgressColumn, console: Console | None = None, auto_refresh: bool = True, refresh_per_second: float = 10, speed_estimate_period: float = 30.0, transient: bool = False, redirect_stdout: bool = True, redirect_stderr: bool = True, get_time: Callable[[], float] | None = None, disable: bool = False, expand: bool = False){:python}`

Bases: `Progress`

Custom Progress bar class, allowing us to have two progress
bars with different columns / layouts.

#### `get_renderables(){:python}`

Get a number of renderables for the progress display.

### _`class{:python}`_`nf_core.pipelines.download.DownloadWorkflow(pipeline=None, revision=None, outdir=None, compress_type=None, force=False, platform=False, download_configuration=None, additional_tags=None, container_system=None, container_library=None, container_cache_utilisation=None, container_cache_index=None, parallel_downloads=4){:python}`

Bases: `object`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

- **Parameters:**
  - **pipeline** (_str_) – A nf-core pipeline name.
  - **revision** (_List_ \*\[\*_str_ _]_) – The workflow revision(s) to download, like 1.0 or dev . Defaults to None.
  - **outdir** (_str_) – Path to the local download directory. Defaults to None.
  - **compress_type** (_str_) – Type of compression for the downloaded files. Defaults to None.
  - **force** (_bool_) – Flag to force download even if files already exist (overwrite existing files). Defaults to False.
  - **platform** (_bool_) – Flag to customize the download for Seqera Platform (convert to git bare repo). Defaults to False.
  - **download_configuration** (_str_) – Download the configuration files from nf-core/configs. Defaults to None.
  - **tag** (_List_ \*\[\*_str_ _]_) – Specify additional tags to add to the downloaded pipeline. Defaults to None.
  - **container_system** (_str_) – The container system to use (e.g., “singularity”). Defaults to None.
  - **container_library** (_List_ \*\[\*_str_ _]_) – The container libraries (registries) to use. Defaults to None.
  - **container_cache_utilisation** (_str_) – If a local or remote cache of already existing container images should be considered. Defaults to None.
  - **container_cache_index** (_str_) – An index for the remote container cache. Defaults to None.
  - **parallel_downloads** (_int_) – The number of parallel downloads to use. Defaults to 4.

#### `compress_download() → None{:python}`

Take the downloaded files and make a compressed .tar.gz archive.

#### `download_configs(){:python}`

Downloads the centralised config profiles from nf-core/configs to `self.outdir`.

#### `download_wf_files(revision, wf_sha, download_url){:python}`

Downloads workflow files from GitHub to the `self.outdir`.

#### `download_workflow(){:python}`

Starts a nf-core workflow download.

#### `download_workflow_platform(location=None){:python}`

Create a bare-cloned git repository of the workflow, so it can be launched with tw launch as <file:/> pipeline

#### `download_workflow_static(){:python}`

Downloads a nf-core workflow from GitHub to the local file system in a self-contained manner.

#### `find_container_images(workflow_directory: str) → None{:python}`

Find container image names for workflow.

Starts by using nextflow config to pull out any process.container
declarations. This works for DSL1. It should return a simple string with resolved logic,
but not always, e.g. not for differentialabundance 1.2.0

Second, we look for DSL2 containers. These can’t be found with
nextflow config at the time of writing, so we scrape the pipeline files.
This returns raw matches that will likely need to be cleaned.

#### `gather_registries(workflow_directory: str) → None{:python}`

Fetch the registries from the pipeline config and CLI arguments and store them in a set.
This is needed to symlink downloaded container images so Nextflow will find them.

#### `get_revision_hash(){:python}`

Find specified revision / branch hash

#### `get_singularity_images(current_revision: str = '') → None{:python}`

Loop through container names and download Singularity images

#### `prioritize_direct_download(container_list){:python}`

Helper function that takes a list of container images (URLs and Docker URIs),
eliminates all Docker URIs for which also a URL is contained and returns the
cleaned and also deduplicated list.

Conceptually, this works like so:

Everything after the last Slash should be identical, e.g. “scanpy:1.7.2–pyhdfd78af_0” in
\[‘https://depot.galaxyproject.org/singularity/scanpy:1.7.2–pyhdfd78af\_0’, ‘biocontainers/scanpy:1.7.2–pyhdfd78af_0’]

re.sub(‘.\*/(.\*)’,’1’,c) will drop everything up to the last slash from c (container_id)

d.get(k:=re.sub(‘.\*/(.\*)’,’1’,c),’’) assigns the truncated string to k (key) and gets the
corresponding value from the dict if present or else defaults to “”.

If the regex pattern matches, the original container_id will be assigned to the dict with the k key.
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

#### `prompt_revision() → None{:python}`

Prompt for pipeline revision / branch
Prompt user for revision tag if ‘–revision’ was not set
If –platform is specified, allow to select multiple revisions
Also the static download allows for multiple revisions, but
we do not prompt this option interactively.

#### `prompt_singularity_cachedir_creation(){:python}`

Prompt about using $NXF_SINGULARITY_CACHEDIR if not already set

#### `prompt_singularity_cachedir_remote(){:python}`

Prompt about the index of a remote $NXF_SINGULARITY_CACHEDIR

#### `prompt_singularity_cachedir_utilization(){:python}`

Ask if we should _only_ use $NXF_SINGULARITY_CACHEDIR without copying into target

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

#### `singularity_copy_cache_image(container: str, out_path: str, cache_path: str | None) → None{:python}`

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

#### `singularity_download_image(container: str, out_path: str, cache_path: str | None, progress:{:python}`[DownloadProgress](#nf_core.pipelines.download.DownloadProgress)) → None

Download a singularity image from the web.

Use native Python to download the file.

- **Parameters:**
  - **container** (_str_) – A pipeline’s container name. Usually it is of similar format
    to `https://depot.galaxyproject.org/singularity/name:version`
  - **out_path** (_str_) – The final target output path
  - **cache_path** (_str_ _,_ _None_) – The NXF_SINGULARITY_CACHEDIR path if set, None if not
  - **progress** (_Progress_) – Rich progress bar instance to add tasks to.

#### `singularity_image_filenames(container: str) → Tuple[str, str | None]{:python}`

Check Singularity cache for image, copy to destination folder if found.

- **Parameters:**
  **container** (_str_) – A pipeline’s container name. Can be direct download URL
  or a Docker Hub repository ID.
- **Returns:**
  Returns a tuple of (out_path, cache_path).
  : out_path is the final target output path. it may point to the NXF_SINGULARITY_CACHEDIR, if cache utilisation was set to ‘amend’.
  If cache utilisation was set to ‘copy’, it will point to the target folder, a subdirectory of the output directory. In the latter case,
  cache_path may either be None (image is not yet cached locally) or point to the image in the NXF_SINGULARITY_CACHEDIR, so it will not be
  downloaded from the web again, but directly copied from there. See get_singularity_images() for implementation.
- **Return type:**
  tuple (str, str)

#### `singularity_pull_image(container: str, out_path: str, cache_path: str | None, library: List[str], progress:{:python}`[DownloadProgress](#nf_core.pipelines.download.DownloadProgress)) → None

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

- **Parameters:**
  - **container** (_str_) – A pipeline’s container name. Usually it is of similar format
    to `nfcore/name:version`.
  - **library** (_list_ _of_ _str_) – A list of libraries to try for pulling the image.
- **Raises:**
  **Various exceptions possible from subprocess execution** **of** **Singularity.** –

#### `symlink_singularity_images(image_out_path: str) → None{:python}`

Create a symlink for each registry in the registry set that points to the image.
We have dropped the explicit registries from the modules in favor of the configurable registries.
Unfortunately, Nextflow still expects the registry to be part of the file name, so a symlink is needed.

The base image, e.g. ./nf-core-gatk-4.4.0.0.img will thus be symlinked as for example ./quay.io-nf-core-gatk-4.4.0.0.img
by prepending all registries in self.registry_set to the image name.

Unfortunately, out output image name may contain a registry definition (Singularity image pulled from depot.galaxyproject.org
or older pipeline version, where the docker registry was part of the image name in the modules). Hence, it must be stripped
before to ensure that it is really the base name.

#### `wf_use_local_configs(revision_dirname){:python}`

Edit the downloaded nextflow.config file to use the local config files

### _`class{:python}`_`nf_core.pipelines.download.WorkflowRepo(remote_url, revision, commit, additional_tags, location=None, hide_progress=False, in_cache=True){:python}`

Bases: `SyncedRepo`

An object to store details about a locally cached workflow repository.

Important Attributes:
: fullname: The full name of the repository, `nf-core/{self.pipelinename}`.
local_repo_dir (str): The local directory, where the workflow is cloned into. Defaults to `$HOME/.cache/nf-core/nf-core/{self.pipeline}`.

#### `__add_additional_tags() → None{:python}`

#### `access(){:python}`

#### `bare_clone(destination){:python}`

#### `checkout(commit){:python}`

Checks out the repository at the requested commit

- **Parameters:**
  **commit** (_str_) – Git SHA of the commit

#### `get_remote_branches(remote_url){:python}`

Get all branches from a remote repository

- **Parameters:**
  **remote_url** (_str_) – The git url to the remote repository
- **Returns:**
  All branches found in the remote
- **Return type:**
  (set\[str])

#### _`property{:python}`_`heads{:python}`

#### `retry_setup_local_repo(skip_confirm=False){:python}`

#### `setup_local_repo(remote, location=None, in_cache=True){:python}`

Sets up the local git repository. If the repository has been cloned previously, it
returns a git.Repo object of that clone. Otherwise it tries to clone the repository from
the provided remote URL and returns a git.Repo of the new clone.

- **Parameters:**
  - **remote** (_str_) – git url of remote
  - **location** (_Path_) – location where the clone should be created/cached.
  - **in_cache** (_bool_ _,_ _optional_) – Whether to clone the repository from the cache. Defaults to False.

Sets self.repo

#### _`property{:python}`_`tags{:python}`

#### `tidy_tags_and_branches(){:python}`

Function to delete all tags and branches that are not of interest to the downloader.
This allows a clutter-free experience in Seqera Platform. The untagged commits are evidently still available.

However, due to local caching, the downloader might also want access to revisions that had been deleted before.
In that case, don’t bother with re-adding the tags and rather download anew from Github.
