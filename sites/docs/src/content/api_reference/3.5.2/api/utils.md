# nf_core.utils

Common utility functions for the nf-core python package.

### _`class{:python}`_`nf_core.utils.GitHubAPISession{:python}`

Bases: `CachedSession`

Class to provide a single session for interacting with the GitHub API for a run.
Inherits the requests_cache.CachedSession and adds additional functionality,
such as automatically setting up GitHub authentication if we can.

#### `get(url, **kwargs){:python}`

Initialise the session if we haven’t already, then call the superclass get method.

#### `lazy_init() → None{:python}`

Initialise the object.

Only do this when it’s actually being used (due to global import)

#### `log_content_headers(request, post_data=None){:python}`

Try to dump everything to the console, useful when things go wrong.

#### `request_retry(url, post_data=None){:python}`

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core pipelines sync code because we get 403 errors: too many simultaneous requests
See <https://github.com/nf-core/tools/issues/911>

#### `safe_get(url){:python}`

Run a GET request, raise a nice exception with lots of logging if it fails.

#### `setup_github_auth(auth=None){:python}`

Try to automatically set up GitHub authentication

### _`class{:python}`_`nf_core.utils.NFCoreTemplateConfig(, org: str | None = None, name: str | None = None, description: str | None = None, author: str | None = None, version: str | None = None, force: bool | None = True, outdir: str | Path | None = None, skip_features: list | None = None, is_nfcore: bool | None = None){:python}`

Bases: `BaseModel`

Template configuration schema

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

#### `author{:python}`_: str | None_

Pipeline author

#### `description{:python}`_: str | None_

Pipeline description

#### `force{:python}`_: bool | None_

Force overwrite of existing files

#### `get(item: str, default: Any = None) → Any{:python}`

#### `is_nfcore{:python}`_: bool | None_

Whether the pipeline is an nf-core pipeline.

#### `model_config{:python}`_= {}_

Configuration for the model, should be a dictionary conforming to \[ConfigDict]\[pydantic.config.ConfigDict].

#### `name{:python}`_: str | None_

Pipeline name

#### `org{:python}`_: str | None_

Organisation name

#### `outdir{:python}`_: str | Path | None_

Output directory

#### _`classmethod{:python}`_`outdir_to_str(v: str | Path | None) → str | None{:python}`

#### `skip_features{:python}`_: list | None_

//nf-co.re/docs/nf-core-tools/pipelines/create for a list of features.

- **Type:**
  Skip features. See https

#### `version{:python}`_: str | None_

Pipeline version

### _`class{:python}`_`nf_core.utils.NFCoreYamlConfig(, repository_type: Literal['pipeline', 'modules'] | None = None, nf_core_version: str | None = None, org_path: str | None = None, lint:{:python}`[`NFCoreYamlLintConfig{:python}`](#nf_core.utils.NFCoreYamlLintConfig)`| None = None, template:{:python}`[`NFCoreTemplateConfig{:python}`](#nf_core.utils.NFCoreTemplateConfig)`| None = None, bump_version: dict[str, bool] | None = None, update: dict[str, str | bool | dict[str, str | dict[str, str | bool]]] | None = None){:python}`

Bases: `BaseModel`

.nf-core.yml configuration file schema

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

#### `bump_version{:python}`_: dict\[str, bool] | None_

//nf-co.re/docs/nf-core-tools/modules/bump-versions for more information.

- **Type:**
  Disable bumping of the version for a module/subworkflow (when repository_type is modules). See https

#### `get(item: str, default: Any = None) → Any{:python}`

#### `lint{:python}`_: [NFCoreYamlLintConfig](#nf_core.utils.NFCoreYamlLintConfig) | None_

//nf-co.re/docs/nf-core-tools/pipelines/lint#linting-config for examples and documentation

- **Type:**
  Pipeline linting configuration, see https

#### `model_config{:python}`_= {}_

Configuration for the model, should be a dictionary conforming to \[ConfigDict]\[pydantic.config.ConfigDict].

#### `model_dump(**kwargs) → dict[str, Any]{:python}`

!!! abstract “Usage Documentation”
: [model_dump](../concepts/serialization#python-mode)

Generate a dictionary representation of the model, optionally specifying which fields to include or exclude.

- **Parameters:**
  - **mode** – The mode in which to_python should run.
    If mode is ‘json’, the output will only contain JSON serializable types.
    If mode is ‘python’, the output may contain non-JSON-serializable Python objects.
  - **include** – A set of fields to include in the output.
  - **exclude** – A set of fields to exclude from the output.
  - **context** – Additional context to pass to the serializer.
  - **by_alias** – Whether to use the field’s alias in the dictionary key if defined.
  - **exclude_unset** – Whether to exclude fields that have not been explicitly set.
  - **exclude_defaults** – Whether to exclude fields that are set to their default value.
  - **exclude_none** – Whether to exclude fields that have a value of None.
  - **exclude_computed_fields** – Whether to exclude computed fields.
    While this can be useful for round-tripping, it is usually recommended to use the dedicated
    round_trip parameter instead.
  - **round_trip** – If True, dumped values should be valid as input for non-idempotent types such as Json\[T].
  - **warnings** – How to handle serialization errors. False/”none” ignores them, True/”warn” logs errors,
    “error” raises a \[PydanticSerializationError]\[pydantic_core.PydanticSerializationError].
  - **fallback** – A function to call when an unknown value is encountered. If not provided,
    a \[PydanticSerializationError]\[pydantic_core.PydanticSerializationError] error is raised.
  - **serialize_as_any** – Whether to serialize fields with duck-typing serialization behavior.
- **Returns:**
  A dictionary representation of the model.

#### `nf_core_version{:python}`_: str | None_

Version of nf-core/tools used to create/update the pipeline

#### `org_path{:python}`_: str | None_

Path to the organisation’s modules repository (used for modules repo_type only)

#### `repository_type{:python}`_: Literal\['pipeline', 'modules'] | None_

Type of repository

#### `template{:python}`_: [NFCoreTemplateConfig](#nf_core.utils.NFCoreTemplateConfig) | None_

Pipeline template configuration

#### `update{:python}`_: dict\[str, str | bool | dict\[str, str | dict\[str, str | bool]]] | None_

//nf-co.re/docs/nf-core-tools/modules/update for more information.

- **Type:**
  Disable updating specific modules/subworkflows (when repository_type is pipeline). See https

### _`class{:python}`_`nf_core.utils.NFCoreYamlLintConfig(, files_unchanged: bool | list[str] | None = None, modules_config: bool | list[str] | None = None, merge_markers: bool | list[str] | None = None, nextflow_config: bool | list[str | dict[str, list[str]]] | None = None, nf_test_content: bool | list[str] | None = None, multiqc_config: bool | list[str] | None = None, files_exist: bool | list[str] | None = None, template_strings: bool | list[str] | None = None, readme: bool | list[str] | None = None, nfcore_components: bool | None = None, actions_nf_test: bool | None = None, actions_awstest: bool | None = None, actions_awsfulltest: bool | None = None, pipeline_todos: bool | None = None, pipeline_if_empty_null: bool | None = None, plugin_includes: bool | None = None, pipeline_name_conventions: bool | None = None, schema_lint: bool | None = None, schema_params: bool | None = None, system_exit: bool | None = None, schema_description: bool | None = None, actions_schema_validation: bool | None = None, modules_json: bool | None = None, modules_structure: bool | None = None, base_config: bool | None = None, nfcore_yml: bool | None = None, version_consistency: bool | None = None, included_configs: bool | None = None, local_component_structure: bool | None = None, rocrate_readme_sync: bool | None = None){:python}`

Bases: `BaseModel`

schema for linting config in .nf-core.yml should cover:

```yaml
files_unchanged:
    - .github/workflows/branch.yml
modules_config: False
modules_config:
        - fastqc
# merge_markers: False
merge_markers:
        - docs/my_pdf.pdf
nextflow_config: False
nextflow_config:
    - manifest.name
    - config_defaults:
        - params.annotation_db
        - params.multiqc_comment_headers
        - params.custom_table_headers
# multiqc_config: False
multiqc_config:
    - report_section_order
    - report_comment
files_exist:
    - .github/CONTRIBUTING.md
    - CITATIONS.md
template_strings: False
template_strings:
        - docs/my_pdf.pdf
nfcore_components: False
# nf_test_content: False
nf_test_content:
    - tests/<test_name>.nf.test
    - tests/nextflow.config
    - nf-test.config
```

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

#### `actions_awsfulltest{:python}`_: bool | None_

Lint all required files to run full tests on AWS

#### `actions_awstest{:python}`_: bool | None_

Lint all required files to run tests on AWS

#### `actions_nf_test{:python}`_: bool | None_

Lint all required files to use GitHub Actions CI

#### `actions_schema_validation{:python}`_: bool | None_

Lint GitHub Action workflow files with schema

#### `base_config{:python}`_: bool | None_

Lint base.config file

#### `files_exist{:python}`_: bool | list\[str] | None_

List of files that can not exist

#### `files_unchanged{:python}`_: bool | list\[str] | None_

List of files that should not be changed

#### `get(item: str, default: Any = None) → Any{:python}`

#### `included_configs{:python}`_: bool | None_

Lint for included configs

#### `local_component_structure{:python}`_: bool | None_

Lint local components use correct structure mirroring remote

#### `merge_markers{:python}`_: bool | list\[str] | None_

List of files that should not contain merge markers

#### `model_config{:python}`_= {}_

Configuration for the model, should be a dictionary conforming to \[ConfigDict]\[pydantic.config.ConfigDict].

#### `modules_config{:python}`_: bool | list\[str] | None_

List of modules that should not be changed

#### `modules_json{:python}`_: bool | None_

Lint modules.json file

#### `modules_structure{:python}`_: bool | None_

Lint modules structure

#### `multiqc_config{:python}`_: bool | list\[str] | None_

List of MultiQC config options that be changed

#### `nextflow_config{:python}`_: bool | list\[str | dict\[str, list\[str]]] | None_

List of Nextflow config files that should not be changed

#### `nf_test_content{:python}`_: bool | list\[str] | None_

List of nf-test content that should not be changed

#### `nfcore_components{:python}`_: bool | None_

Lint all required files to use nf-core modules and subworkflows

#### `nfcore_yml{:python}`_: bool | None_

Lint nf-core.yml

#### `pipeline_if_empty_null{:python}`_: bool | None_

Lint for ifEmpty(null) statements

#### `pipeline_name_conventions{:python}`_: bool | None_

Lint for pipeline name conventions

#### `pipeline_todos{:python}`_: bool | None_

Lint for TODOs statements

#### `plugin_includes{:python}`_: bool | None_

Lint for nextflow plugin

#### `readme{:python}`_: bool | list\[str] | None_

Lint the README.md file

#### `rocrate_readme_sync{:python}`_: bool | None_

Lint for README.md and rocrate.json sync

#### `schema_description{:python}`_: bool | None_

Check that every parameter in the schema has a description.

#### `schema_lint{:python}`_: bool | None_

Lint nextflow_schema.json file

#### `schema_params{:python}`_: bool | None_

Lint schema for all params

#### `system_exit{:python}`_: bool | None_

Lint for System.exit calls in groovy/nextflow code

#### `template_strings{:python}`_: bool | list\[str] | None_

List of files that can contain template strings

#### `version_consistency{:python}`_: bool | None_

Lint for version consistency

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

### _`class{:python}`_`nf_core.utils.SingularityCacheFilePathValidator{:python}`

Bases: `Validator`

Validator for file path specified as –singularity-cache-index argument in nf-core pipelines download

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

#### `validate(value){:python}`

Validate the input.
If invalid, this should raise a `ValidationError`.

- **Parameters:**
  **document** – `Document` instance.

### `nf_core.utils.anaconda_package(dep, dep_channels=None){:python}`

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

### `nf_core.utils.check_nextflow_version(minimal_nf_version: tuple[int, int, int, bool], silent=False) → bool{:python}`

Check the version of Nextflow installed on the system.

- **Parameters:**
  - **minimal_nf_version** (_tuple_ \*\[\*_int_ _,_ _int_ _,_ _int_ _,_ _bool_ _]_) – The minimal version of Nextflow required.
  - **silent** (_bool_) – Whether to log the version or not.
- **Returns:**
  True if the installed version is greater than or equal to minimal_nf_version
- **Return type:**
  bool

### `nf_core.utils.custom_yaml_dumper(){:python}`

Overwrite default PyYAML output to make Prettier YAML linting happy

### `nf_core.utils.determine_base_dir(directory: Path | str = '.') → Path{:python}`

### `nf_core.utils.fetch_remote_version(source_url){:python}`

### `nf_core.utils.fetch_wf_config(wf_path: Path, cache_config: bool = True) → dict{:python}`

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

- **Parameters:**
  - **wf_path** (_str_) – Nextflow workflow file system path.
  - **cache_config** (_bool_) – cache configuration or not (def. True)
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
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(\*\***404\***\*)** –

### `nf_core.utils.get_first_available_path(directory: Path | str, paths: list[str]) → Path | None{:python}`

### `nf_core.utils.get_nf_version() → tuple[int, int, int, bool] | None{:python}`

Get the version of Nextflow installed on the system.

### `nf_core.utils.get_repo_commit(pipeline, commit_id){:python}`

Check if the repo contains the requested commit_id, and expand it to long form if necessary.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **commit_id** – The requested commit ID (SHA). It can be in standard long/short form, or any length.
- **Returns:**
  String or None
- **Return type:**
  commit_id

### `nf_core.utils.get_repo_releases_branches(pipeline, wfs){:python}`

Fetches details of a nf-core workflow to download.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **wfs** – A nf_core.pipelines.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  Array of releases, Array of branches
- **Return type:**
  wf_releases, wf_branches (tuple)
- **Raises:**
  **LockupError**\*\*,\*\* **if the pipeline can not be found.** –

### `nf_core.utils.get_wf_files(wf_path: Path){:python}`

Return a list of all files in a directory (ignores .gitigore files)

### `nf_core.utils.is_file_binary(path){:python}`

Check file path to see if it is a binary file

### `nf_core.utils.is_pipeline_directory(wf_path){:python}`

Checks if the specified directory have the minimum required files
(‘main.nf’, ‘nextflow.config’) for a pipeline directory

- **Parameters:**
  **wf_path** (_str_) – The directory to be inspected
- **Raises:**
  **UserWarning** – If one of the files are missing

### `nf_core.utils.is_relative_to(path1, path2){:python}`

Checks if a path is relative to another.

Should mimic Path.is_relative_to which not available in Python < 3.9

path1 (Path | str): The path that could be a subpath
path2 (Path | str): The path the could be the superpath

### `nf_core.utils.load_tools_config(directory: str | Path = '.') → tuple[Path | None,{:python}`[`NFCoreYamlConfig{:python}`](#nf_core.utils.NFCoreYamlConfig)`| None]{:python}`

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

### `nf_core.utils.plural_es(list_or_int){:python}`

Return a ‘es’ if the input is not one or has not the length of one.

### `nf_core.utils.plural_s(list_or_int){:python}`

Return an s if the input is not one or has not the length of one.

### `nf_core.utils.plural_y(list_or_int){:python}`

Return ‘ies’ if the input is not one or has not the length of one, else ‘y’.

### `nf_core.utils.poll_nfcore_web_api(api_url: str, post_data: dict | None = None) → dict{:python}`

Poll the nf-core website API

Takes argument api_url for URL

Expects API response to be valid JSON and contain a top-level ‘status’ key.

### `nf_core.utils.pretty_nf_version(version: tuple[int, int, int, bool]) → str{:python}`

### `nf_core.utils.prompt_pipeline_release_branch(wf_releases: list[dict[str, Any]], wf_branches: dict[str, Any], multiple: bool = False) → tuple[Any, list[str]]{:python}`

Prompt for pipeline release / branch

- **Parameters:**
  - **wf_releases** (_array_) – Array of repo releases as returned by the GitHub API
  - **wf_branches** (_array_) – Array of repo branches, as returned by the GitHub API
  - **multiple** (_bool_) – Allow selection of multiple releases & branches (for Seqera Platform)
- **Returns:**
  Selected release / branch or False if no releases / branches available
- **Return type:**
  choice (questionary.Choice or bool)

### `nf_core.utils.prompt_remote_pipeline_name(wfs){:python}`

Prompt for the pipeline name with questionary

- **Parameters:**
  **wfs** – A nf_core.pipelines.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  GitHub repo - username/repo
- **Return type:**
  pipeline (str)
- **Raises:**
  **AssertionError**\*\*,\*\* **if pipeline cannot be found** –

### `nf_core.utils.rich_force_colors(){:python}`

Check if any environment variables are set to force Rich to use coloured output

### `nf_core.utils.run_cmd(executable: str, cmd: str) → tuple[bytes, bytes] | None{:python}`

Run a specified command and capture the output. Handle errors nicely.

### `nf_core.utils.set_wd(path: Path) → Generator[None, None, None]{:python}`

Sets the working directory for this context.

- **Parameters:**
  **path** (_Path_) – Path to the working directory to be used inside this context.

### `nf_core.utils.setup_nfcore_cachedir(cache_fn: str | Path) → Path{:python}`

Sets up local caching for caching files between sessions.

### `nf_core.utils.setup_nfcore_dir() → bool{:python}`

Creates a directory for files that need to be kept between sessions

Currently only used for keeping local copies of modules repos

### `nf_core.utils.setup_requests_cachedir() → dict[str, Path | timedelta | str]{:python}`

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .config/nf-core/cache\_\* subdir.

Uses requests_cache monkey patching.
Also returns the config dict so that we can use the same setup with a Session.

### `nf_core.utils.sort_dictionary(d: dict) → dict{:python}`

Sorts a nested dictionary recursively

### `nf_core.utils.strip_ansi_codes(string, replace_with=''){:python}`

Strip ANSI colouring codes from a string to return plain text.

From Stack Overflow: <https://stackoverflow.com/a/14693789/713980>

### `nf_core.utils.unquote(s: str) → str{:python}`

Remove paired quotes (single or double) from start and end of string.

Uses ast.literal_eval to safely parse Python string literals, preserving
the original string if it’s not a valid literal.

- **Parameters:**
  **s** – String potentially containing quotes
- **Returns:**
  String with outer quotes removed if present, otherwise original string

### `nf_core.utils.validate_file_md5(file_name, expected_md5hex){:python}`

Validates the md5 checksum of a file on disk.

- **Parameters:**
  - **file_name** (_str_) – Path to a local file.
  - **expected** (_str_) – The expected md5sum.
- **Raises:**
  **IOError**\*\*,\*\* **if the md5sum does not match the remote sum.** –

### `nf_core.utils.wait_cli_function(poll_func: Callable[[], bool], refresh_per_second: int = 20) → None{:python}`

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

- **Parameters:**
  - **poll_func** (_function_) – Function to call
  - **refresh_per_second** (_int_) – Refresh this many times per second. Default: 20.
- **Returns:**
  None. Just sits in an infinite loop until the function returns True.
