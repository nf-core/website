<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.utils`

Common utility functions for the nf-core python package.

---

<a href="../../../../../../tools/nf_core/utils.py#L49"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L72"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `rich_force_colors`

```python
rich_force_colors()
```

Check if any environment variables are set to force Rich to use coloured output

---

<a href="../../../../../../tools/nf_core/utils.py#L81"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `github_api_auto_auth`

```python
github_api_auto_auth()
```

---

<a href="../../../../../../tools/nf_core/utils.py#L184"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetch_wf_config`

```python
fetch_wf_config(wf_path)
```

Uses Nextflow to retrieve the the configuration variables from a Nextflow workflow.

**Args:**

- <b>`wf_path`</b> (str): Nextflow workflow file system path.

**Returns:**

- <b>`dict`</b>: Workflow configuration settings.

---

<a href="../../../../../../tools/nf_core/utils.py#L264"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `nextflow_cmd`

```python
nextflow_cmd(cmd)
```

Run a Nextflow command and capture the output. Handle errors nicely

---

<a href="../../../../../../tools/nf_core/utils.py#L278"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `setup_requests_cachedir`

```python
setup_requests_cachedir()
```

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user's home directory under a .nfcore_cache subdir.

---

<a href="../../../../../../tools/nf_core/utils.py#L295"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L319"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `poll_nfcore_web_api`

```python
poll_nfcore_web_api(api_url, post_data=None)
```

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level 'status' key.

---

<a href="../../../../../../tools/nf_core/utils.py#L359"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L411"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_anaconda_licence`

```python
parse_anaconda_licence(anaconda_response, version=None)
```

Given a response from the anaconda API using anaconda_package, parse the software licences.

Returns: Set of licence types

---

<a href="../../../../../../tools/nf_core/utils.py#L445"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L472"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_biocontainer_tag`

```python
get_biocontainer_tag(package, version)
```

Given a bioconda package and version, look for a container at quay.io and returns the tag of the most recent image that matches the package version Sends a HTTP GET request to the quay.io API.

**Args:**

- <b>`package`</b> (str): A bioconda package name.
- <b>`version`</b> (str): Version of the bioconda package

**Raises:**
A LookupError, if the connection fails or times out or gives an unexpected status code A ValueError, if the package name can not be found (404)

---

<a href="../../../../../../tools/nf_core/utils.py#L519"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `custom_yaml_dumper`

```python
custom_yaml_dumper()
```

Overwrite default PyYAML output to make Prettier YAML linting happy

---

<a href="../../../../../../tools/nf_core/utils.py#L92"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/utils.py#L110"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(wf_path)
```

Initialise pipeline object

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
