# nf_core.utils

Common utility functions for the nf-core python package.

### _`class{:python}`_`nf_core.utils.GitHub_API_Session{:python}`

Bases: `CachedSession`

Class to provide a single session for interacting with the GitHub API for a run.
Inherits the requests_cache.CachedSession and adds additional functionality,
such as automatically setting up GitHub authentication if we can.

#### `get(url, **kwargs){:python}`

Initialise the session if we haven’t already, then call the superclass get method.

#### `lazy_init(){:python}`

Initialise the object.

Only do this when it’s actually being used (due to global import)

#### `log_content_headers(request){:python}`

Try to dump everything to the console, useful when things go wrong.

#### `request_retry(url, post_data=None){:python}`

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core sync code because we get 403 errors: too many simultaneous requests
See <https://github.com/nf-core/tools/issues/911>

#### `safe_get(url){:python}`

Run a GET request, raise a nice exception with lots of logging if it fails.

#### `setup_github_auth(auth=None){:python}`

Try to automatically set up GitHub authentication

### _`class{:python}`_`nf_core.utils.Pipeline(wf_path){:python}`

Bases: `object`

Object to hold information about a local pipeline.

- **Parameters:**
  **path** (_str_) – The path to the nf-core pipeline directory.

#### `conda_config{:python}`

The parsed conda configuration file content (`environment.yml`).

- **Type:**
  dict

#### `conda_package_info{:python}`

The conda package(s) information, based on the API requests to Anaconda cloud.

- **Type:**
  dict

#### `nf_config{:python}`

The Nextflow pipeline configuration file content.

- **Type:**
  dict

#### `files{:python}`

A list of files found during the linting process.

- **Type:**
  list

#### `git_sha{:python}`

The git sha for the repo commit / current GitHub pull-request ($GITHUB_PR_COMMIT)

- **Type:**
  str

#### `minNextflowVersion{:python}`

The minimum required Nextflow version to run the pipeline.

- **Type:**
  str

#### `wf_path{:python}`

Path to the pipeline directory.

- **Type:**
  str

#### `pipeline_name{:python}`

The pipeline name, without the nf-core tag, for example hlatyping.

- **Type:**
  str

#### `schema_obj{:python}`

A `PipelineSchema` object

- **Type:**
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

### `nf_core.utils.anaconda_package(dep, dep_channels=['conda-forge', 'bioconda', 'defaults']){:python}`

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

- **Parameters:**
  - **dep** (_str_) – A conda package name.
  - **dep_channels** (_list_) – list of conda channels to use
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(\*\***404\***\*)** –

### `nf_core.utils.check_if_outdated(current_version=None, remote_version=None, source_url='https://nf-co.re/tools_version'){:python}`

Check if the current version of nf-core is outdated

### `nf_core.utils.custom_yaml_dumper(){:python}`

Overwrite default PyYAML output to make Prettier YAML linting happy

### `nf_core.utils.fetch_wf_config(wf_path, cache_config=True){:python}`

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

- **Parameters:**
  - **wf_path** (_str_) – Nextflow workflow file system path.
  - **cache_config** (_bool_) – cache configuration or not (def. True)
- **Returns:**
  Workflow configuration settings.
- **Return type:**
  dict

### `nf_core.utils.get_biocontainer_tag(package, version){:python}`

Given a bioconda package and version, looks for Docker and Singularity containers
using the biocontaineres API, e.g.:
<https://api.biocontainers.pro/ga4gh/trs/v2/tools>/{tool}/versions/{tool}-{version}
Returns the most recent container versions by default.
:param package: A bioconda package name.
:type package: str
:param version: Version of the bioconda package
:type version: str

- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(\*\***404\***\*)** –

### `nf_core.utils.get_repo_releases_branches(pipeline, wfs){:python}`

Fetches details of a nf-core workflow to download.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **wfs** – A nf_core.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  Array of releases, Array of branches
- **Return type:**
  wf_releases, wf_branches (tuple)
- **Raises:**
  **LockupError**\*\*,\*\* **if the pipeline can not be found.** –

### `nf_core.utils.is_file_binary(path){:python}`

Check file path to see if it is a binary file

### `nf_core.utils.is_pipeline_directory(wf_path){:python}`

Checks if the specified directory have the minimum required files
(‘main.nf’, ‘nextflow.config’) for a pipeline directory

- **Parameters:**
  **wf_path** (_str_) – The directory to be inspected
- **Raises:**
  **UserWarning** – If one of the files are missing

### `nf_core.utils.load_tools_config(dir='.'){:python}`

Parse the nf-core.yml configuration file

Look for a file called either .nf-core.yml or .nf-core.yaml

Also looks for the deprecated file .nf-core-lint.yml/yaml and issues
a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn’t be loaded

### `nf_core.utils.nextflow_cmd(cmd){:python}`

Run a Nextflow command and capture the output. Handle errors nicely

### `nf_core.utils.parse_anaconda_licence(anaconda_response, version=None){:python}`

Given a response from the anaconda API using anaconda_package, parse the software licences.

Returns: Set of licence types

### `nf_core.utils.pip_package(dep){:python}`

Query PyPI package information.

Sends a HTTP GET request to the PyPI remote API.

- **Parameters:**
  **dep** (_str_) – A PyPI package name.
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** –

### `nf_core.utils.poll_nfcore_web_api(api_url, post_data=None){:python}`

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level ‘status’ key.

### `nf_core.utils.prompt_pipeline_release_branch(wf_releases, wf_branches){:python}`

Prompt for pipeline release / branch

- **Parameters:**
  - **wf_releases** (_array_) – Array of repo releases as returned by the GitHub API
  - **wf_branches** (_array_) – Array of repo branches, as returned by the GitHub API
- **Returns:**
  Selected release / branch name
- **Return type:**
  choice (str)

### `nf_core.utils.prompt_remote_pipeline_name(wfs){:python}`

Prompt for the pipeline name with questionary

- **Parameters:**
  **wfs** – A nf_core.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  GitHub repo - username/repo
- **Return type:**
  pipeline (str)
- **Raises:**
  **AssertionError**\*\*,\*\* **if pipeline cannot be found** –

### `nf_core.utils.rich_force_colors(){:python}`

Check if any environment variables are set to force Rich to use coloured output

### `nf_core.utils.setup_requests_cachedir(){:python}`

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .config/nf-core/cache\_\* subdir.

Uses requests_cache monkey patching.
Also returns the config dict so that we can use the same setup with a Session.

### `nf_core.utils.sort_dictionary(d){:python}`

Sorts a nested dictionary recursively

### `nf_core.utils.wait_cli_function(poll_func, poll_every=20){:python}`

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

- **Parameters:**
  - **poll_func** (_function_) – Function to call
  - **poll_every** (_int_) – How many tenths of a second to wait between function calls. Default: 20.
- **Returns:**
  None. Just sits in an infite loop until the function returns True.
