# nf\_core.utils

Common utility functions for the nf-core python package.

### *`class{:python}`*`nf_core.utils.GitHubAPISession{:python}`

Bases: `CachedSession`

Class to provide a single session for interacting with the GitHub API for a run.
Inherits the requests\_cache.CachedSession and adds additional functionality,
such as automatically setting up GitHub authentication if we can.

#### `get(url, **kwargs){:python}`

Initialise the session if we haven’t already, then call the superclass get method.

#### `lazy_init(){:python}`

Initialise the object.

Only do this when it’s actually being used (due to global import)

#### `log_content_headers(request, post_data=None){:python}`

Try to dump everything to the console, useful when things go wrong.

#### `request_retry(url, post_data=None){:python}`

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core sync code because we get 403 errors: too many simultaneous requests
See <https://github.com/nf-core/tools/issues/911>

#### `safe_get(url){:python}`

Run a GET request, raise a nice exception with lots of logging if it fails.

#### `setup_github_auth(auth=None){:python}`

Try to automatically set up GitHub authentication

### *`class{:python}`*`nf_core.utils.Pipeline(wf_path){:python}`

Bases: `object`

Object to hold information about a local pipeline.

* **Parameters:**
  **path** (*str*) – The path to the nf-core pipeline directory.

#### `conda_config{:python}`

The parsed conda configuration file content (`environment.yml`).

* **Type:**
  dict

#### `conda_package_info{:python}`

The conda package(s) information, based on the API requests to Anaconda cloud.

* **Type:**
  dict

#### `nf_config{:python}`

The Nextflow pipeline configuration file content.

* **Type:**
  dict

#### `files{:python}`

A list of files found during the linting process.

* **Type:**
  list

#### `git_sha{:python}`

The git sha for the repo commit / current GitHub pull-request ($GITHUB\_PR\_COMMIT)

* **Type:**
  str

#### `minNextflowVersion{:python}`

The minimum required Nextflow version to run the pipeline.

* **Type:**
  str

#### `wf_path{:python}`

Path to the pipeline directory.

* **Type:**
  str

#### `pipeline_name{:python}`

The pipeline name, without the nf-core tag, for example hlatyping.

* **Type:**
  str

#### `schema_obj{:python}`

A `PipelineSchema` object

* **Type:**
  obj

#### `_fp(fn){:python}`

Convenience function to get full path to a file in the pipeline

#### `_list_files(){:python}`

Get a list of all files in the pipeline

#### `_load(){:python}`

Run core load functions

#### `_load_conda_environment(){:python}`

Try to load the pipeline environment.yml file, if it exists

#### `_load_pipeline_config(){:python}`

Get the nextflow config for this pipeline

Once loaded, set a few convienence reference class attributes

### *`class{:python}`*`nf_core.utils.SingularityCacheFilePathValidator{:python}`

Bases: `Validator`

Validator for file path specified as –singularity-cache-index argument in nf-core download

#### `_abc_impl{:python}`*= <\_abc.\_abc\_data object>*

#### `validate(value){:python}`

Validate the input.
If invalid, this should raise a `ValidationError`.

* **Parameters:**
  **document** – `Document` instance.

### `nf_core.utils.anaconda_package(dep, dep_channels=None){:python}`

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

* **Parameters:**
  * **dep** (*str*) – A conda package name.
  * **dep\_channels** (*list*) – list of conda channels to use
* **Raises:**
  * **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  * **A ValueError**\*\*,\*\* **if the package name can not be found** **(****404****)** –

### `nf_core.utils.check_if_outdated(current_version=None, remote_version=None, source_url='https://nf-co.re/tools_version'){:python}`

Check if the current version of nf-core is outdated

### `nf_core.utils.custom_yaml_dumper(){:python}`

Overwrite default PyYAML output to make Prettier YAML linting happy

### `nf_core.utils.determine_base_dir(directory='.'){:python}`

### `nf_core.utils.fetch_remote_version(source_url){:python}`

### `nf_core.utils.fetch_wf_config(wf_path, cache_config=True){:python}`

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

* **Parameters:**
  * **wf\_path** (*str*) – Nextflow workflow file system path.
  * **cache\_config** (*bool*) – cache configuration or not (def. True)
* **Returns:**
  Workflow configuration settings.
* **Return type:**
  dict

### `nf_core.utils.file_md5(fname){:python}`

Calculates the md5sum for a file on the disk.

* **Parameters:**
  **fname** (*str*) – Path to a local file.

### `nf_core.utils.get_biocontainer_tag(package, version){:python}`

Given a bioconda package and version, looks for Docker and Singularity containers
using the biocontaineres API, e.g.:
<https://api.biocontainers.pro/ga4gh/trs/v2/tools>/{tool}/versions/{tool}-{version}
Returns the most recent container versions by default.
:param package: A bioconda package name.
:type package: str
:param version: Version of the bioconda package
:type version: str

* **Raises:**
  * **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  * **A ValueError**\*\*,\*\* **if the package name can not be found** **(****404****)** –

### `nf_core.utils.get_first_available_path(directory, paths){:python}`

### `nf_core.utils.get_repo_releases_branches(pipeline, wfs){:python}`

Fetches details of a nf-core workflow to download.

* **Parameters:**
  * **pipeline** (*str*) – GitHub repo username/repo
  * **wfs** – A nf\_core.list.Workflows() object, where get\_remote\_workflows() has been called.
* **Returns:**
  Array of releases, Array of branches
* **Return type:**
  wf\_releases, wf\_branches (tuple)
* **Raises:**
  **LockupError**\*\*,\*\* **if the pipeline can not be found.** –

### `nf_core.utils.is_file_binary(path){:python}`

Check file path to see if it is a binary file

### `nf_core.utils.is_pipeline_directory(wf_path){:python}`

Checks if the specified directory have the minimum required files
(‘main.nf’, ‘nextflow.config’) for a pipeline directory

* **Parameters:**
  **wf\_path** (*str*) – The directory to be inspected
* **Raises:**
  **UserWarning** – If one of the files are missing

### `nf_core.utils.is_relative_to(path1, path2){:python}`

Checks if a path is relative to another.

Should mimic Path.is\_relative\_to which not available in Python < 3.9

path1 (Path | str): The path that could be a subpath
path2 (Path | str): The path the could be the superpath

### `nf_core.utils.load_tools_config(directory: str | Path = '.'){:python}`

Parse the nf-core.yml configuration file

Look for a file called either .nf-core.yml or .nf-core.yaml

Also looks for the deprecated file .nf-core-lint.yml/yaml and issues
a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn’t be loaded

### `nf_core.utils.nested_delitem(d, keys){:python}`

Deletes a key from a nested dictionary

* **Parameters:**
  * **d** (*dict*) – the nested dictionary to traverse
  * **keys** (*list* \*\[\**Any* *]*) – A list of keys to iteratively traverse, deleting the final one

### `nf_core.utils.nested_setitem(d, keys, value){:python}`

Sets the value in a nested dict using a list of keys to traverse

* **Parameters:**
  * **d** (*dict*) – the nested dictionary to traverse
  * **keys** (*list* \*\[\**Any* *]*) – A list of keys to iteratively traverse
  * **value** (*Any*) – The value to be set for the last key in the chain

### `nf_core.utils.parse_anaconda_licence(anaconda_response, version=None){:python}`

Given a response from the anaconda API using anaconda\_package, parse the software licences.

Returns: Set of licence types

### `nf_core.utils.pip_package(dep){:python}`

Query PyPI package information.

Sends a HTTP GET request to the PyPI remote API.

* **Parameters:**
  **dep** (*str*) – A PyPI package name.
* **Raises:**
  * **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** –
  * **A ValueError**\*\*,\*\* **if the package name can not be found** –

### `nf_core.utils.plural_es(list_or_int){:python}`

Return a ‘es’ if the input is not one or has not the length of one.

### `nf_core.utils.plural_s(list_or_int){:python}`

Return an s if the input is not one or has not the length of one.

### `nf_core.utils.plural_y(list_or_int){:python}`

Return ‘ies’ if the input is not one or has not the length of one, else ‘y’.

### `nf_core.utils.poll_nfcore_web_api(api_url, post_data=None){:python}`

Poll the nf-core website API

Takes argument api\_url for URL

Expects API reponse to be valid JSON and contain a top-level ‘status’ key.

### `nf_core.utils.prompt_pipeline_release_branch(wf_releases, wf_branches, multiple=False){:python}`

Prompt for pipeline release / branch

* **Parameters:**
  * **wf\_releases** (*array*) – Array of repo releases as returned by the GitHub API
  * **wf\_branches** (*array*) – Array of repo branches, as returned by the GitHub API
  * **multiple** (*bool*) – Allow selection of multiple releases & branches (for Tower)
* **Returns:**
  Selected release / branch name
* **Return type:**
  choice (str)

### `nf_core.utils.prompt_remote_pipeline_name(wfs){:python}`

Prompt for the pipeline name with questionary

* **Parameters:**
  **wfs** – A nf\_core.list.Workflows() object, where get\_remote\_workflows() has been called.
* **Returns:**
  GitHub repo - username/repo
* **Return type:**
  pipeline (str)
* **Raises:**
  **AssertionError**\*\*,\*\* **if pipeline cannot be found** –

### `nf_core.utils.rich_force_colors(){:python}`

Check if any environment variables are set to force Rich to use coloured output

### `nf_core.utils.run_cmd(executable: str, cmd: str){:python}`

Run a specified command and capture the output. Handle errors nicely.

### `nf_core.utils.set_wd(path: Path){:python}`

Sets the working directory for this context.

* **Parameters:**
  **path** (*Path*) – Path to the working directory to be used inside this context.

### `nf_core.utils.setup_nfcore_cachedir(cache_fn: str | Path){:python}`

Sets up local caching for caching files between sessions.

### `nf_core.utils.setup_nfcore_dir(){:python}`

Creates a directory for files that need to be kept between sessions

Currently only used for keeping local copies of modules repos

### `nf_core.utils.setup_requests_cachedir(){:python}`

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .config/nf-core/cache\_\* subdir.

Uses requests\_cache monkey patching.
Also returns the config dict so that we can use the same setup with a Session.

### `nf_core.utils.sort_dictionary(d){:python}`

Sorts a nested dictionary recursively

### `nf_core.utils.strip_ansi_codes(string, replace_with=''){:python}`

Strip ANSI colouring codes from a string to return plain text.

From Stack Overflow: <https://stackoverflow.com/a/14693789/713980>

### `nf_core.utils.validate_file_md5(file_name, expected_md5hex){:python}`

Validates the md5 checksum of a file on disk.

* **Parameters:**
  * **file\_name** (*str*) – Path to a local file.
  * **expected** (*str*) – The expected md5sum.
* **Raises:**
  **IOError**\*\*,\*\* **if the md5sum does not match the remote sum.** –

### `nf_core.utils.wait_cli_function(poll_func, refresh_per_second=20){:python}`

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

* **Parameters:**
  * **poll\_func** (*function*) – Function to call
  * **refresh\_per\_second** (*int*) – Refresh this many times per second. Default: 20.
* **Returns:**
  None. Just sits in an infite loop until the function returns True.
