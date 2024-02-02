# nf_core.download

Downloads a nf-core pipeline to the local file system.

### _exception_ nf_core.download.ContainerError(container, registry, address, absolute_URI, out_path, singularity_command, error_msg)

Bases: `Exception`

A class of errors related to pulling containers with Singularity/Apptainer

#### _exception_ ImageExists(error_log)

Bases: `FileExistsError`

Image already exists in cache/output directory.

#### _exception_ ImageNotFound(error_log)

Bases: `FileNotFoundError`

The image can not be found in the registry

#### _exception_ InvalidTag(error_log)

Bases: `AttributeError`

Image and registry are valid, but the (version) tag is not

#### _exception_ OtherError(error_log)

Bases: `RuntimeError`

Undefined error with the container

#### _exception_ RegistryNotFound(error_log)

Bases: `ConnectionRefusedError`

The specified registry does not resolve to a valid IP address

### _exception_ nf_core.download.DownloadError

Bases: `RuntimeError`

A custom exception that is raised when nf-core download encounters a problem that we already took into consideration.
In this case, we do not want to print the traceback, but give the user some concise, helpful feedback instead.

### _class_ nf_core.download.DownloadProgress(\*columns: str | ProgressColumn, console: Console | None = None, auto_refresh: bool = True, refresh_per_second: float = 10, speed_estimate_period: float = 30.0, transient: bool = False, redirect_stdout: bool = True, redirect_stderr: bool = True, get_time: Callable[[], float] | None = None, disable: bool = False, expand: bool = False)

Bases: `Progress`

Custom Progress bar class, allowing us to have two progress
bars with different columns / layouts.

#### get_renderables()

Get a number of renderables for the progress display.

### _class_ nf_core.download.DownloadWorkflow(pipeline=None, revision=None, outdir=None, compress_type=None, force=False, tower=False, download_configuration=None, container_system=None, container_library=None, container_cache_utilisation=None, container_cache_index=None, parallel_downloads=4)

Bases: `object`

Downloads a nf-core workflow from GitHub to the local file system.

Can also download its Singularity container image if required.

- **Parameters:**
  - **pipeline** (_str_) – A nf-core pipeline name.
  - **revision** (_List**[**str\*\*]_) – The workflow revision to download, like 1.0. Defaults to None.
  - **container** (_bool_) – Flag, if the Singularity container should be downloaded as well. Defaults to False.
  - **tower** (_bool_) – Flag, to customize the download for Nextflow Tower (convert to git bare repo). Defaults to False.
  - **outdir** (_str_) – Path to the local download directory. Defaults to None.

#### compress_download()

Take the downloaded files and make a compressed .tar.gz archive.

#### download_configs()

Downloads the centralised config profiles from nf-core/configs to `self.outdir`.

#### download_wf_files(revision, wf_sha, download_url)

Downloads workflow files from GitHub to the `self.outdir`.

#### download_workflow()

Starts a nf-core workflow download.

#### download_workflow_static()

Downloads a nf-core workflow from GitHub to the local file system in a self-contained manner.

#### download_workflow_tower(location=None)

Create a bare-cloned git repository of the workflow, so it can be launched with tw launch as [file:/](file:/) pipeline

#### find_container_images(workflow_directory)

Find container image names for workflow.

Starts by using nextflow config to pull out any process.container
declarations. This works for DSL1. It should return a simple string with resolved logic,
but not always, e.g. not for differentialabundance 1.2.0

Second, we look for DSL2 containers. These can’t be found with
nextflow config at the time of writing, so we scrape the pipeline files.
This returns raw matches that will likely need to be cleaned.

#### get_revision_hash()

Find specified revision / branch hash

#### get_singularity_images(current_revision='')

Loop through container names and download Singularity images

#### prioritize_direct_download(container_list)

Helper function that takes a list of container images (URLs and Docker URIs),
eliminates all Docker URIs for which also a URL is contained and returns the
cleaned and also deduplicated list.

Conceptually, this works like so:

Everything after the last Slash should be identical, e.g. “scanpy:1.7.2–pyhdfd78af_0” in
[‘https://depot.galaxyproject.org/singularity/scanpy:1.7.2–pyhdfd78af_0’, ‘biocontainers/scanpy:1.7.2–pyhdfd78af_0’]

re.sub(‘.\*/(.\*)’,’1’,c) will drop everything up to the last slash from c (container_id)

d.get(k:=re.sub(‘.\*/(.\*)’,’1’,c),’’) assigns the truncated string to k (key) and gets the
corresponding value from the dict if present or else defaults to “”.

If the regex pattern matches, the original container_id will be assigned to the dict with the k key.
r”^$|(?!^http)” matches an empty string (we didn’t have it in the dict yet and want to keep it in either case) or
any string that does not start with http. Because if our current dict value already starts with http,
we want to keep it and not replace with with whatever we have now (which might be the Docker URI).

A regex that matches http, r”^$|^http” could thus be used to prioritize the Docker URIs over http Downloads

#### prompt_compression_type()

Ask user if we should compress the downloaded files

#### prompt_config_inclusion()

Prompt for inclusion of institutional configurations

#### prompt_container_download()

Prompt whether to download container images or not

#### prompt_pipeline_name()

Prompt for the pipeline name if not set with a flag

#### prompt_revision()

Prompt for pipeline revision / branch
Prompt user for revision tag if ‘–revision’ was not set
If –tower is specified, allow to select multiple revisions
Also the static download allows for multiple revisions, but
we do not prompt this option interactively.

#### prompt_singularity_cachedir_creation()

Prompt about using $NXF_SINGULARITY_CACHEDIR if not already set

#### prompt_singularity_cachedir_remote()

Prompt about the index of a remote $NXF_SINGULARITY_CACHEDIR

#### prompt_singularity_cachedir_utilization()

Ask if we should _only_ use $NXF_SINGULARITY_CACHEDIR without copying into target

#### read_remote_containers()

Reads the file specified as index for the remote Singularity cache dir

#### rectify_raw_container_matches(raw_findings)

Helper function to rectify the raw extracted container matches into fully qualified container names.
If multiple containers are found, any prefixed with http for direct download is prioritized

Example syntax:

Early DSL2:
: if (workflow.containerEngine == ‘singularity’ && !params.singularity_pull_docker_container) {
: container “[https://depot.galaxyproject.org/singularity/fastqc:0.11.9–0](https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0)”
<br/>
} else {
: container “quay.io/biocontainers/fastqc:0.11.9–0”
<br/>
}

Later DSL2:
: container “${ workflow.containerEngine == ‘singularity’ && !task.ext.singularity_pull_docker_container ?
: ‘[https://depot.galaxyproject.org/singularity/fastqc:0.11.9–0](https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0)’ :
‘biocontainers/fastqc:0.11.9–0’ }”

Later DSL2, variable is being used:
: container “${ workflow.containerEngine == ‘singularity’ && !task.ext.singularity_pull_docker_container ?
  : “[https://depot.galaxyproject.org/singularity](https://depot.galaxyproject.org/singularity)/${container_id}” :
“quay.io/biocontainers/${container_id}” }”
<br/>
container_id = ‘mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0’

DSL1 / Special case DSL2:
: container “nfcore/cellranger:6.0.2”

#### singularity_copy_cache_image(container, out_path, cache_path)

Copy Singularity image from NXF_SINGULARITY_CACHEDIR to target folder.

#### singularity_download_image(container, out_path, cache_path, progress)

Download a singularity image from the web.

Use native Python to download the file.

- **Parameters:**
  - **container** (_str_) – A pipeline’s container name. Usually it is of similar format
    to `https://depot.galaxyproject.org/singularity/name:version`
  - **out_path** (_str_) – The final target output path
  - **cache_path** (_str\*\*,_ _None_) – The NXF_SINGULARITY_CACHEDIR path if set, None if not
  - **progress** (_Progress_) – Rich progress bar instance to add tasks to.

#### singularity_image_filenames(container)

Check Singularity cache for image, copy to destination folder if found.

- **Parameters:**
  **container** (_str_) – A pipeline’s container name. Can be direct download URL
  or a Docker Hub repository ID.
- **Returns:**
  Returns True if we have the image in the target location.
  : Returns a download path if not.
- **Return type:**
  results (bool, str)

#### singularity_pull_image(container, out_path, cache_path, library, progress)

Pull a singularity image using `singularity pull`

Attempt to use a local installation of singularity to pull the image.

- **Parameters:**
  - **container** (_str_) – A pipeline’s container name. Usually it is of similar format
    to `nfcore/name:version`.
  - **library** (_list_ _of_ _str_) – A list of libraries to try for pulling the image.
- **Raises:**
  **Various exceptions possible from subprocess execution** **of** **Singularity.** –

#### wf_use_local_configs(revision_dirname)

Edit the downloaded nextflow.config file to use the local config files

### _class_ nf_core.download.WorkflowRepo(remote_url, revision, commit, location=None, hide_progress=False, in_cache=True)

Bases: `SyncedRepo`

An object to store details about a locally cached workflow repository.

Important Attributes:
: fullname: The full name of the repository, `nf-core/{self.pipelinename}`.
local_repo_dir (str): The local directory, where the workflow is cloned into. Defaults to `$HOME/.cache/nf-core/nf-core/{self.pipeline}`.

#### access()

#### bare_clone(destination)

#### checkout(commit)

Checks out the repository at the requested commit

- **Parameters:**
  **commit** (_str_) – Git SHA of the commit

#### get_remote_branches(remote_url)

Get all branches from a remote repository

- **Parameters:**
  **remote_url** (_str_) – The git url to the remote repository
- **Returns:**
  All branches found in the remote
- **Return type:**
  (set[str])

#### retry_setup_local_repo(skip_confirm=False)

#### setup_local_repo(remote, location=None, in_cache=True)

Sets up the local git repository. If the repository has been cloned previously, it
returns a git.Repo object of that clone. Otherwise it tries to clone the repository from
the provided remote URL and returns a git.Repo of the new clone.

- **Parameters:**
  - **remote** (_str_) – git url of remote
  - **location** (_Path_) – location where the clone should be created/cached.
  - **in_cache** (_bool\*\*,_ _optional_) – Whether to clone the repository from the cache. Defaults to False.

Sets self.repo

#### tidy_tags_and_branches()

Function to delete all tags and branches that are not of interest to the downloader.
This allows a clutter-free experience in Tower. The untagged commits are evidently still available.

However, due to local caching, the downloader might also want access to revisions that had been deleted before.
In that case, don’t bother with re-adding the tags and rather download anew from Github.
