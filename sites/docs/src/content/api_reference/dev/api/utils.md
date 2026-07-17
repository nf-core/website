# nf\_core.utils

Common utility functions for the nf-core python package.

### _`class{:python}`_`nf_core.utils.ContainerRegistryUrls(*values){:python}`

Bases: `Enum`

#### `GALAXY_SINGULARITY{:python}`_= 'depot.galaxyproject.org/singularity'_

#### `SEQERA_DOCKER{:python}`_= 'community.wave.seqera.io/library'_

#### `SEQERA_SINGULARITY{:python}`_= 'community-cr-prod.seqera.io/docker/registry/v2'_

### _`class{:python}`_`nf_core.utils.Pipeline(wf_path: Path){:python}`

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

The git sha for the repo commit / current GitHub pull-request ($GITHUB\_PR\_COMMIT)

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

#### `_fp(fn: str | Path) → Path{:python}`

Convenience function to get full path to a file in the pipeline

#### `_load() → bool{:python}`

Run core load functions

#### `_load_conda_environment() → bool{:python}`

Try to load the pipeline environment.yml file, if it exists

#### `conda_config{:python}`_: dict_

#### `conda_package_info{:python}`_: dict_

#### `files{:python}`_: list\[Path]_

#### `git_sha{:python}`_: str | None_

#### `list_files() → list[Path]{:python}`

Get a list of all files in the pipeline

#### `load_pipeline_config() → bool{:python}`

Get the nextflow config for this pipeline

Once loaded, set a few convenience reference class attributes

#### `minNextflowVersion{:python}`_: str | None_

#### `nf_config{:python}`_: dict_

#### `pipeline_name{:python}`_: str | None_

#### `pipeline_prefix{:python}`_: str | None_

#### `repo{:python}`_: git.Repo | None_

#### `schema_obj{:python}`_: [PipelineSchema](pipelines/schema#nf_core.pipelines.schema.PipelineSchema) | None_

### `nf_core.utils._load_cached_remote_version() → str | None{:python}`

### `nf_core.utils._load_version_cache() → dict{:python}`

### `nf_core.utils._nfcore_question_style(){:python}`

### `nf_core.utils._spawn_remote_version_refresh(source_url: str) → None{:python}`

Kick off a detached background process to refresh the remote version cache.

### `nf_core.utils.anaconda_package(dep, dep_channels=None){:python}`

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

- **Parameters:**
  - **dep** (_str_) – A conda package name.
  - **dep\_channels** (_list_) – list of conda channels to use
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(****404****)** –

### `nf_core.utils.check_if_outdated(current_version=None, remote_version=None, source_url='https://nf-co.re/tools_version'){:python}`

Check if the current version of nf-core is outdated

### `nf_core.utils.check_nextflow_version(minimal_nf_version: tuple[int, int, int, bool], silent=False) → bool{:python}`

Check the version of Nextflow installed on the system.

- **Parameters:**
  - **minimal\_nf\_version** (_tuple_ \*\[\*_int_ _,_ _int_ _,_ _int_ _,_ _bool_ _]_) – The minimal version of Nextflow required.
  - **silent** (_bool_) – Whether to log the version or not.
- **Returns:**
  True if the installed version is greater than or equal to minimal\_nf\_version
- **Return type:**
  bool

### `nf_core.utils.custom_yaml_dumper(){:python}`

Overwrite default PyYAML output to make Prettier YAML linting happy

### `nf_core.utils.determine_base_dir(directory: Path | str = '.') → Path{:python}`

### `nf_core.utils.fetch_wf_config(wf_path: Path, cache_config: bool = True) → dict{:python}`

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

- **Parameters:**
  - **wf\_path** (_str_) – Nextflow workflow file system path.
  - **cache\_config** (_bool_) – cache configuration or not (def. True)
- **Returns:**
  Workflow configuration settings.
- **Return type:**
  dict

### `nf_core.utils.file_md5(fname){:python}`

Calculates the md5sum for a file on the disk.

- **Parameters:**
  **fname** (_str_) – Path to a local file.

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
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(****404****)** –

### `nf_core.utils.get_first_available_path(directory: Path | str, paths: list[str]) → Path | None{:python}`

### `nf_core.utils.get_nf_version() → tuple[int, int, int, bool] | None{:python}`

Get the version of Nextflow installed on the system. Cached for the lifetime of the process.

### `nf_core.utils.get_repo_commit(pipeline, commit_id){:python}`

Check if the repo contains the requested commit\_id, and expand it to long form if necessary.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **commit\_id** – The requested commit ID (SHA). It can be in standard long/short form, or any length.
- **Returns:**
  String or None
- **Return type:**
  commit\_id

### `nf_core.utils.get_repo_releases_branches(pipeline, wfs){:python}`

Fetches details of a nf-core workflow to download.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **wfs** – A nf\_core.pipelines.list.Workflows() object, where get\_remote\_workflows() has been called.
- **Returns:**
  Array of releases, Array of branches
- **Return type:**
  wf\_releases, wf\_branches (tuple)
- **Raises:**
  **LockupError**\*\*,\*\* **if the pipeline can not be found.** –

### `nf_core.utils.get_wf_files(wf_path: Path){:python}`

Return a list of all files in a directory (ignores .gitigore files)

### `nf_core.utils.is_file_binary(path){:python}`

Check file path to see if it is a binary file

### `nf_core.utils.is_interactive() → bool{:python}`

Check if the current session is interactive (has a TTY on stdin, stdout, and stderr).

### `nf_core.utils.is_pipeline_directory(wf_path){:python}`

Checks if the specified directory have the minimum required files
(‘main.nf’, ‘nextflow.config’) for a pipeline directory

- **Parameters:**
  **wf\_path** (_str_) – The directory to be inspected
- **Raises:**
  **UserWarning** – If one of the files are missing

### `nf_core.utils.is_relative_to(path1, path2){:python}`

Checks if a path is relative to another.

Should mimic Path.is\_relative\_to which not available in Python < 3.9

path1 (Path | str): The path that could be a subpath
path2 (Path | str): The path the could be the superpath

### `nf_core.utils.lazy_attrs(module_globals: dict, mapping: dict[str, str | Callable]){:python}`

Build module-level `__getattr__` and `__dir__` for lazy attributes (PEP 562).

Each attribute maps either to the dotted path of the module to import it
from, or to a zero-argument callable that builds it. Resolved attributes
are cached in the module’s globals, so the hook only fires once per name.

Usage:

```default
__getattr__, __dir__ = lazy_attrs(globals(), {"PipelineCreateApp": "nf_core.pipelines.create"})
```

### `nf_core.utils.load_tools_config(directory: str | Path = '.') → tuple[Path | None,{:python}`[`NFCoreYamlConfig{:python}`](nf_core_yml#nf_core.utils.NFCoreYamlConfig)`| None]{:python}`

Parse the nf-core.yml configuration file

Look for a file called either .nf-core.yml or .nf-core.yaml

Also looks for the deprecated file .nf-core-lint.yml/yaml and issues
a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn’t be loaded

### `nf_core.utils.nested_delitem(d, keys){:python}`

Deletes a key from a nested dictionary

- **Parameters:**
  - **d** (_dict_) – the nested dictionary to traverse
  - **keys** (_list_ \*\[\*_Any_ _]_) – A list of keys to iteratively traverse, deleting the final one

### `nf_core.utils.nested_setitem(d, keys, value){:python}`

Sets the value in a nested dict using a list of keys to traverse

- **Parameters:**
  - **d** (_dict_) – the nested dictionary to traverse
  - **keys** (_list_ \*\[\*_Any_ _]_) – A list of keys to iteratively traverse
  - **value** (_Any_) – The value to be set for the last key in the chain

### `nf_core.utils.parse_anaconda_licence(anaconda_response, version=None){:python}`

Given a response from the anaconda API using anaconda\_package, parse the software licences.

Returns: Set of licence types

### `nf_core.utils.pip_package(dep){:python}`

Query PyPI package information.

Sends a HTTP GET request to the PyPI remote API.

- **Parameters:**
  **dep** (_str_) – A PyPI package name.
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** –

### `nf_core.utils.plural_es(list_or_int){:python}`

Return a ‘es’ if the input is not one or has not the length of one.

### `nf_core.utils.plural_s(list_or_int){:python}`

Return an s if the input is not one or has not the length of one.

### `nf_core.utils.plural_y(list_or_int){:python}`

Return ‘ies’ if the input is not one or has not the length of one, else ‘y’.

### `nf_core.utils.poll_nfcore_web_api(api_url: str, post_data: dict | None = None) → dict{:python}`

Poll the nf-core website API

Takes argument api\_url for URL

Expects API response to be valid JSON and contain a top-level ‘status’ key.

### `nf_core.utils.pretty_nf_version(version: tuple[int, int, int, bool]) → str{:python}`

### `nf_core.utils.prompt_pipeline_release_branch(wf_releases: list[dict[str, Any]], wf_branches: dict[str, Any], multiple: bool = False) → tuple[Any, list[str]]{:python}`

Prompt for pipeline release / branch

- **Parameters:**
  - **wf\_releases** (_array_) – Array of repo releases as returned by the GitHub API
  - **wf\_branches** (_array_) – Array of repo branches, as returned by the GitHub API
  - **multiple** (_bool_) – Allow selection of multiple releases & branches (for Seqera Platform)
- **Returns:**
  Selected release / branch or False if no releases / branches available
- **Return type:**
  choice (questionary.Choice or bool)

### `nf_core.utils.prompt_remote_pipeline_name(wfs){:python}`

Prompt for the pipeline name with questionary

- **Parameters:**
  **wfs** – A nf\_core.pipelines.list.Workflows() object, where get\_remote\_workflows() has been called.
- **Returns:**
  GitHub repo - username/repo
- **Return type:**
  pipeline (str)
- **Raises:**
  **AssertionError**\*\*,\*\* **if pipeline cannot be found** –

### `nf_core.utils.read_module_name(main_nf: Path) → str | None{:python}`

Return the process name declared in a Nextflow `main.nf` file, or `None`.

### `nf_core.utils.rich_force_colors(){:python}`

Check if any environment variables are set to force Rich to use coloured output

### `nf_core.utils.run_cmd(executable: str, cmd: str) → tuple[bytes, bytes] | None{:python}`

Run a specified command and capture the output. Handle errors nicely.

### `nf_core.utils.set_wd(path: Path) → Generator[None, None, None]{:python}`

Sets the working directory for this context.

- **Parameters:**
  **path** (_Path_) – Path to the working directory to be used inside this context.

### `nf_core.utils.set_wd_tempdir(base_dir: Path | None = None) → Generator[Path, None, None]{:python}`

Context manager to create a tempdir and change into it, ensuring its removal and a return to
the original working directory on exit (including exceptions).

- **Parameters:**
  **base\_dir** – Directory in which to create the tempdir. Defaults to the system temp location.

### `nf_core.utils.setup_nfcore_cachedir(cache_fn: str | Path | None = None) → Path{:python}`

Sets up local caching for caching files between sessions.

Creates and returns a subdirectory of the nf-core cache directory,
or the cache directory itself if no subdirectory is given.

### `nf_core.utils.setup_requests_cachedir() → dict[str, Path | timedelta | str]{:python}`

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .config/nf-core/cache\_\* subdir.

Uses requests\_cache monkey patching.
Also returns the config dict so that we can use the same setup with a Session.

### `nf_core.utils.sort_dictionary(d: dict) → dict{:python}`

Sorts a nested dictionary recursively

### `nf_core.utils.strip_ansi_codes(string, replace_with=''){:python}`

Strip ANSI colouring codes from a string to return plain text.

From Stack Overflow: <https://stackoverflow.com/a/14693789/713980>

### `nf_core.utils.unquote(s: str) → str{:python}`

Remove paired quotes (single or double) from start and end of string.

Uses ast.literal\_eval to safely parse Python string literals, preserving
the original string if it’s not a valid literal.

Special handling for ruamel.yaml DoubleQuotedScalarString to preserve
strings that should not be converted to numbers (e.g., “123” stays as string).

- **Parameters:**
  **s** – String potentially containing quotes
- **Returns:**
  String with outer quotes removed if present, otherwise original string

### `nf_core.utils.validate_file_md5(file_name, expected_md5hex){:python}`

Validates the md5 checksum of a file on disk.

- **Parameters:**
  - **file\_name** (_str_) – Path to a local file.
  - **expected** (_str_) – The expected md5sum.
- **Raises:**
  **IOError**\*\*,\*\* **if the md5sum does not match the remote sum.** –

### `nf_core.utils.wait_cli_function(poll_func: Callable[[], bool], refresh_per_second: int = 20) → None{:python}`

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

- **Parameters:**
  - **poll\_func** (_function_) – Function to call
  - **refresh\_per\_second** (_int_) – Refresh this many times per second. Default: 20.
- **Returns:**
  None. Just sits in an infinite loop until the function returns True.
