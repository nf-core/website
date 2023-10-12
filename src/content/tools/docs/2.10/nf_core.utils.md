<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.utils`

Common utility functions for the nf-core python package.

## **Global Variables**

- **NFCORE_CACHE_DIR**
- **NFCORE_DIR**
- **gh_api**
- **CONFIG_PATHS**
- **DEPRECATED_CONFIG_PATHS**

---

<a href="../../../../../../tools/nf_core/utils.py#L62"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetch_remote_version`

```python
fetch_remote_version(source_url)
```

---

<a href="../../../../../../tools/nf_core/utils.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_if_outdated`

```python
check_if_outdated(
    current_version=None,
    remote_version=None,
    source_url='https://nf-co.re/tools_version'
)
```

Check if the current version of nf-core is outdated

---

<a href="../../../../../../tools/nf_core/utils.py#L97"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `rich_force_colors`

```python
rich_force_colors()
```

Check if any environment variables are set to force Rich to use coloured output

---

<a href="../../../../../../tools/nf_core/utils.py#L199"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `is_pipeline_directory`

```python
is_pipeline_directory(wf_path)
```

Checks if the specified directory have the minimum required files ('main.nf', 'nextflow.config') for a pipeline directory

**Args:**

- <b>`wf_path`</b> (str): The directory to be inspected

**Raises:**

- <b>`UserWarning`</b>: If one of the files are missing

---

<a href="../../../../../../tools/nf_core/utils.py#L216"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetch_wf_config`

```python
fetch_wf_config(wf_path, cache_config=True)
```

Uses Nextflow to retrieve the the configuration variables from a Nextflow workflow.

**Args:**

- <b>`wf_path`</b> (str): Nextflow workflow file system path.
- <b>`cache_config`</b> (bool): cache configuration or not (def. True)

**Returns:**

- <b>`dict`</b>: Workflow configuration settings.

---

<a href="../../../../../../tools/nf_core/utils.py#L306"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `nextflow_cmd`

```python
nextflow_cmd(cmd)
```

Run a Nextflow command and capture the output. Handle errors nicely

---

<a href="../../../../../../tools/nf_core/utils.py#L320"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_nfcore_dir`

```python
setup_nfcore_dir()
```

Creates a directory for files that need to be kept between sessions

Currently only used for keeping local copies of modules repos

---

<a href="../../../../../../tools/nf_core/utils.py#L329"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_requests_cachedir`

```python
setup_requests_cachedir()
```

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user's home directory under a .config/nf-core/cache\_\* subdir.

Uses requests_cache monkey patching. Also returns the config dict so that we can use the same setup with a Session.

---

<a href="../../../../../../tools/nf_core/utils.py#L358"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `wait_cli_function`

```python
wait_cli_function(poll_func, refresh_per_second=20)
```

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

**Arguments:**

- <b>`poll_func`</b> (function): Function to call
- <b>`refresh_per_second`</b> (int): Refresh this many times per second. Default: 20.

**Returns:**
None. Just sits in an infite loop until the function returns True.

---

<a href="../../../../../../tools/nf_core/utils.py#L382"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `poll_nfcore_web_api`

```python
poll_nfcore_web_api(api_url, post_data=None)
```

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level 'status' key.

---

<a href="../../../../../../tools/nf_core/utils.py#L572"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `anaconda_package`

```python
anaconda_package(dep, dep_channels=None)
```

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

**Args:**

- <b>`dep`</b> (str): A conda package name.
- <b>`dep_channels`</b> (list): list of conda channels to use

**Raises:**
A LookupError, if the connection fails or times out or gives an unexpected status code A ValueError, if the package name can not be found (404)

---

<a href="../../../../../../tools/nf_core/utils.py#L626"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_anaconda_licence`

```python
parse_anaconda_licence(anaconda_response, version=None)
```

Given a response from the anaconda API using anaconda_package, parse the software licences.

Returns: Set of licence types

---

<a href="../../../../../../tools/nf_core/utils.py#L660"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `pip_package`

```python
pip_package(dep)
```

Query PyPI package information.

Sends a HTTP GET request to the PyPI remote API.

**Args:**

- <b>`dep`</b> (str): A PyPI package name.

**Raises:**
A LookupError, if the connection fails or times out A ValueError, if the package name can not be found

---

<a href="../../../../../../tools/nf_core/utils.py#L686"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_biocontainer_tag`

```python
get_biocontainer_tag(package, version)
```

Given a bioconda package and version, looks for Docker and Singularity containers using the biocontaineres API, e.g.: https://api.biocontainers.pro/ga4gh/trs/v2/tools/{tool}/versions/{tool}-{version} Returns the most recent container versions by default.

**Args:**

- <b>`package`</b> (str): A bioconda package name.
- <b>`version`</b> (str): Version of the bioconda package

**Raises:**
A LookupError, if the connection fails or times out or gives an unexpected status code A ValueError, if the package name can not be found (404)

---

<a href="../../../../../../tools/nf_core/utils.py#L754"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `custom_yaml_dumper`

```python
custom_yaml_dumper()
```

Overwrite default PyYAML output to make Prettier YAML linting happy

---

<a href="../../../../../../tools/nf_core/utils.py#L786"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `is_file_binary`

```python
is_file_binary(path)
```

Check file path to see if it is a binary file

---

<a href="../../../../../../tools/nf_core/utils.py#L802"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_remote_pipeline_name`

```python
prompt_remote_pipeline_name(wfs)
```

Prompt for the pipeline name with questionary

**Args:**

- <b>`wfs`</b>: A nf_core.list.Workflows() object, where get_remote_workflows() has been called.

**Returns:**

- <b>`pipeline`</b> (str): GitHub repo - username/repo

**Raises:**
AssertionError, if pipeline cannot be found

---

<a href="../../../../../../tools/nf_core/utils.py#L840"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_pipeline_release_branch`

```python
prompt_pipeline_release_branch(wf_releases, wf_branches, multiple=False)
```

Prompt for pipeline release / branch

**Args:**

- <b>`wf_releases`</b> (array): Array of repo releases as returned by the GitHub API
- <b>`wf_branches`</b> (array): Array of repo branches, as returned by the GitHub API
- <b>`multiple`</b> (bool): Allow selection of multiple releases & branches (for Tower)

**Returns:**

- <b>`choice`</b> (str): Selected release / branch name

---

<a href="../../../../../../tools/nf_core/utils.py#L901"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_repo_releases_branches`

```python
get_repo_releases_branches(pipeline, wfs)
```

Fetches details of a nf-core workflow to download.

**Args:**

- <b>`pipeline`</b> (str): GitHub repo username/repo
- <b>`wfs`</b>: A nf_core.list.Workflows() object, where get_remote_workflows() has been called.

**Returns:**

- <b>`wf_releases, wf_branches (tuple)`</b>: Array of releases, Array of branches

**Raises:**
LockupError, if the pipeline can not be found.

---

<a href="../../../../../../tools/nf_core/utils.py#L978"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `load_tools_config`

```python
load_tools_config(directory='.')
```

Parse the nf-core.yml configuration file

Look for a file called either `.nf-core.yml` or `.nf-core.yaml`

Also looks for the deprecated file `.nf-core-lint.yml/yaml` and issues a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn't be loaded

---

<a href="../../../../../../tools/nf_core/utils.py#L1014"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `determine_base_dir`

```python
determine_base_dir(directory='.')
```

---

<a href="../../../../../../tools/nf_core/utils.py#L1024"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_first_available_path`

```python
get_first_available_path(directory, paths)
```

---

<a href="../../../../../../tools/nf_core/utils.py#L1031"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `sort_dictionary`

```python
sort_dictionary(d)
```

Sorts a nested dictionary recursively

---

<a href="../../../../../../tools/nf_core/utils.py#L1042"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `plural_s`

```python
plural_s(list_or_int)
```

Return an s if the input is not one or has not the length of one.

---

<a href="../../../../../../tools/nf_core/utils.py#L1048"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `plural_y`

```python
plural_y(list_or_int)
```

Return 'ies' if the input is not one or has not the length of one, else 'y'.

---

<a href="../../../../../../tools/nf_core/utils.py#L1054"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `plural_es`

```python
plural_es(list_or_int)
```

Return a 'es' if the input is not one or has not the length of one.

---

<a href="../../../../../../tools/nf_core/utils.py#L1065"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `strip_ansi_codes`

```python
strip_ansi_codes(string, replace_with='')
```

Strip ANSI colouring codes from a string to return plain text.

From Stack Overflow: https://stackoverflow.com/a/14693789/713980

---

<a href="../../../../../../tools/nf_core/utils.py#L1073"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `is_relative_to`

```python
is_relative_to(path1, path2)
```

Checks if a path is relative to another.

Should mimic Path.is_relative_to which not available in Python < 3.9

path1 (Path | str): The path that could be a subpath path2 (Path | str): The path the could be the superpath

---

<a href="../../../../../../tools/nf_core/utils.py#L1085"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `file_md5`

```python
file_md5(fname)
```

Calculates the md5sum for a file on the disk.

**Args:**

- <b>`fname`</b> (str): Path to a local file.

---

<a href="../../../../../../tools/nf_core/utils.py#L1101"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `validate_file_md5`

```python
validate_file_md5(file_name, expected_md5hex)
```

Validates the md5 checksum of a file on disk.

**Args:**

- <b>`file_name`</b> (str): Path to a local file.
- <b>`expected`</b> (str): The expected md5sum.

**Raises:**
IOError, if the md5sum does not match the remote sum.

---

<a href="../../../../../../tools/nf_core/utils.py#L106"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Pipeline`

Object to hold information about a local pipeline.

**Args:**

- <b>`path`</b> (str): The path to the nf-core pipeline directory.

**Attributes:**

- <b>`conda_config`</b> (dict): The parsed conda configuration file content (`environment.yml`).
- <b>`conda_package_info`</b> (dict): The conda package(s) information, based on the API requests to Anaconda cloud.
- <b>`nf_config`</b> (dict): The Nextflow pipeline configuration file content.
- <b>`files`</b> (list): A list of files found during the linting process.
- <b>`git_sha`</b> (str): The git sha for the repo commit / current GitHub pull-request (`$GITHUB_PR_COMMIT`)
- <b>`minNextflowVersion`</b> (str): The minimum required Nextflow version to run the pipeline.
- <b>`wf_path`</b> (str): Path to the pipeline directory.
- <b>`pipeline_name`</b> (str): The pipeline name, without the `nf-core` tag, for example `hlatyping`.
- <b>`schema_obj`</b> (obj): A :class:`PipelineSchema` object

<a href="../../../../../../tools/nf_core/utils.py#L124"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(wf_path)
```

Initialise pipeline object

---

<a href="../../../../../../tools/nf_core/utils.py#L421"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `GitHub_API_Session`

Class to provide a single session for interacting with the GitHub API for a run. Inherits the requests_cache.CachedSession and adds additional functionality, such as automatically setting up GitHub authentication if we can.

<a href="../../../../../../tools/nf_core/utils.py#L428"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

---

<a href="../../../../../../tools/nf_core/utils.py#L518"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get`

```python
get(url, **kwargs)
```

Initialise the session if we haven't already, then call the superclass get method.

---

<a href="../../../../../../tools/nf_core/utils.py#L434"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lazy_init`

```python
lazy_init()
```

Initialise the object.

Only do this when it's actually being used (due to global import)

---

<a href="../../../../../../tools/nf_core/utils.py#L486"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `log_content_headers`

```python
log_content_headers(request, post_data=None)
```

Try to dump everything to the console, useful when things go wrong.

---

<a href="../../../../../../tools/nf_core/utils.py#L526"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `request_retry`

```python
request_retry(url, post_data=None)
```

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core sync code because we get 403 errors: too many simultaneous requests See https://github.com/nf-core/tools/issues/911

---

<a href="../../../../../../tools/nf_core/utils.py#L506"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `safe_get`

```python
safe_get(url)
```

Run a GET request, raise a nice exception with lots of logging if it fails.

---

<a href="../../../../../../tools/nf_core/utils.py#L446"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_github_auth`

```python
setup_github_auth(auth=None)
```

Try to automatically set up GitHub authentication

---

<a href="../../../../../../tools/nf_core/utils.py#L884"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SingularityCacheFilePathValidator`

Validator for file path specified as --singularity-cache-index argument in nf-core download

---

<a href="../../../../../../tools/nf_core/utils.py#L889"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate(value)
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
