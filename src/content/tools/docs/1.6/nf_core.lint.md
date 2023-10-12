<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint`

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to the nf-core community guidelines.

---

<a href="../../../../../../tools/nf_core/lint.py#L28"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `run_linting`

```python
run_linting(pipeline_dir, release_mode=False)
```

Runs all nf-core linting checks on a given Nextflow pipeline project in either `release` mode or `normal` mode (default). Returns an object of type :class:`PipelineLint` after finished.

**Args:**

- <b>`pipeline_dir`</b> (str): The path to the Nextflow pipeline root directory
- <b>`release_mode`</b> (bool): Set this to `True`, if the linting should be run in the `release` mode.
- <b>`See `</b>: class:`PipelineLint` for more information.

**Returns:**

- <b>`An object of type `</b>: class:`PipelineLint` that contains all the linting results.

---

<a href="../../../../../../tools/nf_core/lint.py#L67"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineLint`

Object to hold linting information and results. All objects attributes are set, after the :func:`PipelineLint.lint_pipeline` function was called.

**Args:**

- <b>`path`</b> (str): The path to the nf-core pipeline directory.

**Attributes:**

- <b>`conda_config`</b> (dict): The parsed conda configuration file content (`environment.yml`).
- <b>`conda_package_info`</b> (dict): The conda package(s) information, based on the API requests to Anaconda cloud.
- <b>`config`</b> (dict): The Nextflow pipeline configuration file content.
- <b>`dockerfile`</b> (list): A list of lines (str) from the parsed Dockerfile.
- <b>`failed`</b> (list): A list of tuples of the form: `(<error no>, <reason>)`
- <b>`files`</b> (list): A list of files found during the linting process.
- <b>`minNextflowVersion`</b> (str): The minimum required Nextflow version to run the pipeline.
- <b>`passed`</b> (list): A list of tuples of the form: `(<passed no>, <reason>)`
- <b>`path`</b> (str): Path to the pipeline directory.
- <b>`pipeline_name`</b> (str): The pipeline name, without the `nf-core` tag, for example `hlatyping`.
- <b>`release_mode`</b> (bool): `True`, if you the to linting was run in release mode, `False` else.
- <b>`warned`</b> (list): A list of tuples of the form: `(<warned no>, <reason>)`

**Attribute specifications**

Some of the more complex attributes of a PipelineLint object.

- `conda_config`:

```

    # Example       {          'name': 'nf-core-hlatyping',          'channels': ['bioconda', 'conda-forge'],          'dependencies': ['optitype=1.3.2', 'yara=0.9.6']       }

* `conda_package_info`:
```

    # See https://api.anaconda.org/package/bioconda/bioconda-utils as an example.       {          <package>: <API JSON repsonse object>       }

```
* `config`: Produced by calling Nextflow with :code:`nextflow config -flat <workflow dir>`. Here is an example from the `nf-core/hlatyping <https://github.com/nf-core/hlatyping>`_ pipeline:
```

         process.container = 'nfcore/hlatyping:1.1.1'          params.help = false          params.outdir = './results'          params.bam = false          params.singleEnd = false          params.seqtype = 'dna'          params.solver = 'glpk'          params.igenomes_base = './iGenomes'          params.clusterOptions = false          ...

<a href="../../../../../../tools/nf_core/lint.py#L122"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(path)
```

Initialise linting object

---

<a href="../../../../../../tools/nf_core/lint.py#L658"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_anaconda_package`

```python
check_anaconda_package(dep)
```

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

**Args:**

- <b>`dep`</b> (str): A conda package name.

**Raises:**
A ValueError, if the package name can not be resolved.

---

<a href="../../../../../../tools/nf_core/lint.py#L439"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_ci_config`

```python
check_ci_config()
```

Checks that the Travis or Circle CI YAML config is valid.

Makes sure that `nf-core lint` runs in travis tests and that tests run with the required nextflow version.

---

<a href="../../../../../../tools/nf_core/lint.py#L740"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_conda_dockerfile`

```python
check_conda_dockerfile()
```

Checks the Docker build file.

Checks that: _ a name is given and is consistent with the pipeline name _ dependency versions are pinned \* dependency versions are the latest available

---

<a href="../../../../../../tools/nf_core/lint.py#L585"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_conda_env_yaml`

```python
check_conda_env_yaml()
```

Checks that the conda environment file is valid.

Checks that: _ a name is given and is consistent with the pipeline name _ check that dependency versions are pinned \* dependency versions are the latest available

---

<a href="../../../../../../tools/nf_core/lint.py#L266"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_docker`

```python
check_docker()
```

Checks that Dockerfile contains the string `FROM`.

---

<a href="../../../../../../tools/nf_core/lint.py#L188"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_files_exist`

```python
check_files_exist()
```

Checks a given pipeline directory for required files.

Iterates through the pipeline's directory content and checkmarks files for presence. Files that **must** be present:

```

     'nextflow.config',      'Dockerfile',      ['.travis.yml', '.circle.yml'],      ['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling      'README.md',      'CHANGELOG.md',      'docs/README.md',      'docs/output.md',      'docs/usage.md'

Files that *should* be present:
```

     'main.nf',      'environment.yml',      'conf/base.config'

````


**Raises:**
  An AssertionError if neither `nextflow.config` or `main.nf` found.

---

<a href="../../../../../../tools/nf_core/lint.py#L280"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_licence`

```python
check_licence()
````

Checks licence file is MIT.

Currently the checkpoints are: _ licence file must be long enough (4 or more lines) _ licence contains the string _without restriction_ \* licence doesn't have any placeholder variables

---

<a href="../../../../../../tools/nf_core/lint.py#L322"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_nextflow_config`

```python
check_nextflow_config()
```

Checks a given pipeline for required config variables.

Uses `nextflow config -flat` to parse pipeline `nextflow.config` and print all config variables. NB: Does NOT parse contents of main.nf / nextflow script

---

<a href="../../../../../../tools/nf_core/lint.py#L711"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_pip_package`

```python
check_pip_package(dep)
```

Query PyPi package information.

Sends a HTTP GET request to the PyPi remote API.

**Args:**

- <b>`dep`</b> (str): A PyPi package name.

**Raises:**
A ValueError, if the package name can not be resolved or the connection timed out.

---

<a href="../../../../../../tools/nf_core/lint.py#L765"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_pipeline_todos`

```python
check_pipeline_todos()
```

Go through all template files looking for the string 'TODO nf-core:'

---

<a href="../../../../../../tools/nf_core/lint.py#L505"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_readme`

```python
check_readme()
```

Checks the repository README file for errors.

Currently just checks the badges at the top of the README.

---

<a href="../../../../../../tools/nf_core/lint.py#L537"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_version_consistency`

```python
check_version_consistency()
```

Checks container tags versions.

Runs on `process.container`, `process.container` and `$TRAVIS_TAG` (each only if set).

Checks that: _ the container has a tag _ the version numbers are numeric \* the version numbers are the same as one-another

---

<a href="../../../../../../tools/nf_core/lint.py#L137"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_pipeline`

```python
lint_pipeline(release_mode=False)
```

Main linting function.

Takes the pipeline directory as the primary input and iterates through the different linting checks in order. Collects any warnings or errors and returns summary at completion. Raises an exception if there is a critical error that makes the rest of the tests pointless (eg. no pipeline script). Results from this function are printed by the main script.

**Args:**

- <b>`release_mode`</b> (boolean): Activates the release mode, which checks for consistent version tags of containers. Default is `False`.

**Returns:**
dict: Summary of test result messages structured as follows:

```

    {          'pass': [              ( test-id (int), message (string) ),              ( test-id (int), message (string) )          ],          'warn': [(id, msg)],          'fail': [(id, msg)],     }

```

**Raises:**
If a critical problem is found, an `AssertionError` is raised.

---

<a href="../../../../../../tools/nf_core/lint.py#L788"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_results`

```python
print_results()
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
