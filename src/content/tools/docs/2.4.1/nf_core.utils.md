<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.utils`

Common utility functions for the nf-core python package.

## **Global Variables**

- **NFCORE_CONFIG_DIR**
- **gh_api**

---

<a href="../../../../../../tools/nf_core/utils.py#L58"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L81"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `rich_force_colors`

```python
rich_force_colors()
```

Check if any environment variables are set to force Rich to use coloured output

---

<a href="../../../../../../tools/nf_core/utils.py#L182"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L199"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L286"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `nextflow_cmd`

```python
nextflow_cmd(cmd)
```

Run a Nextflow command and capture the output. Handle errors nicely

---

<a href="../../../../../../tools/nf_core/utils.py#L300"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_requests_cachedir`

```python
setup_requests_cachedir()
```

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user's home directory under a .config/nf-core/cache\_\* subdir.

Uses requests_cache monkey patching. Also returns the config dict so that we can use the same setup with a Session.

---

<a href="../../../../../../tools/nf_core/utils.py#L328"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `wait_cli_function`

```python
wait_cli_function(poll_func, poll_every=20)
```

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

**Arguments:**

- <b>`poll_func`</b> (function): Function to call
- <b>`poll_every`</b> (int): How many tenths of a second to wait between function calls. Default: 20.

**Returns:**
None. Just sits in an infite loop until the function returns True.

---

<a href="../../../../../../tools/nf_core/utils.py#L352"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `poll_nfcore_web_api`

```python
poll_nfcore_web_api(api_url, post_data=None)
```

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level 'status' key.

---

<a href="../../../../../../tools/nf_core/utils.py#L543"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `anaconda_package`

```python
anaconda_package(dep, dep_channels=['conda-forge', 'bioconda', 'defaults'])
```

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

**Args:**

- <b>`dep`</b> (str): A conda package name.
- <b>`dep_channels`</b> (list): list of conda channels to use

**Raises:**
A LookupError, if the connection fails or times out or gives an unexpected status code A ValueError, if the package name can not be found (404)

---

<a href="../../../../../../tools/nf_core/utils.py#L595"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_anaconda_licence`

```python
parse_anaconda_licence(anaconda_response, version=None)
```

Given a response from the anaconda API using anaconda_package, parse the software licences.

Returns: Set of licence types

---

<a href="../../../../../../tools/nf_core/utils.py#L629"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L656"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L708"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `custom_yaml_dumper`

```python
custom_yaml_dumper()
```

Overwrite default PyYAML output to make Prettier YAML linting happy

---

<a href="../../../../../../tools/nf_core/utils.py#L740"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `is_file_binary`

```python
is_file_binary(path)
```

Check file path to see if it is a binary file

---

<a href="../../../../../../tools/nf_core/utils.py#L756"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L795"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_pipeline_release_branch`

```python
prompt_pipeline_release_branch(wf_releases, wf_branches)
```

Prompt for pipeline release / branch

**Args:**

- <b>`wf_releases`</b> (array): Array of repo releases as returned by the GitHub API
- <b>`wf_branches`</b> (array): Array of repo branches, as returned by the GitHub API

**Returns:**

- <b>`choice`</b> (str): Selected release / branch name

---

<a href="../../../../../../tools/nf_core/utils.py#L825"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L901"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `load_tools_config`

```python
load_tools_config(dir='.')
```

Parse the nf-core.yml configuration file

Look for a file called either `.nf-core.yml` or `.nf-core.yaml`

Also looks for the deprecated file `.nf-core-lint.yml/yaml` and issues a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn't be loaded

---

<a href="../../../../../../tools/nf_core/utils.py#L941"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `sort_dictionary`

```python
sort_dictionary(d)
```

Sorts a nested dictionary recursively

---

<a href="../../../../../../tools/nf_core/utils.py#L90"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L108"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(wf_path)
```

Initialise pipeline object

---

<a href="../../../../../../tools/nf_core/utils.py#L392"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `GitHub_API_Session`

Class to provide a single session for interacting with the GitHub API for a run. Inherits the requests_cache.CachedSession and adds additional functionality, such as automatically setting up GitHub authentication if we can.

<a href="../../../../../../tools/nf_core/utils.py#L399"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

---

<a href="../../../../../../tools/nf_core/utils.py#L489"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get`

```python
get(url, **kwargs)
```

Initialise the session if we haven't already, then call the superclass get method.

---

<a href="../../../../../../tools/nf_core/utils.py#L405"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lazy_init`

```python
lazy_init()
```

Initialise the object.

Only do this when it's actually being used (due to global import)

---

<a href="../../../../../../tools/nf_core/utils.py#L457"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `log_content_headers`

```python
log_content_headers(request, post_data=None)
```

Try to dump everything to the console, useful when things go wrong.

---

<a href="../../../../../../tools/nf_core/utils.py#L497"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `request_retry`

```python
request_retry(url, post_data=None)
```

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core sync code because we get 403 errors: too many simultaneous requests See https://github.com/nf-core/tools/issues/911

---

<a href="../../../../../../tools/nf_core/utils.py#L477"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `safe_get`

```python
safe_get(url)
```

Run a GET request, raise a nice exception with lots of logging if it fails.

---

<a href="../../../../../../tools/nf_core/utils.py#L417"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_github_auth`

```python
setup_github_auth(auth=None)
```

Try to automatically set up GitHub authentication

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
