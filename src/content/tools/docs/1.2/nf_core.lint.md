<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint`

Linting code for the nf-core python package.

Tests Nextflow pipelines to check that they adhere to the nf-core community guidelines.

## **Global Variables**

- **cachedir**

---

<a href="../../../../../../tools/nf_core/lint.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `run_linting`

```python
run_linting(pipeline_dir, release)
```

Run all linting tests. Called by main script.

---

<a href="../../../../../../tools/nf_core/lint.py#L65"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineLint`

Object to hold linting info and results

<a href="../../../../../../tools/nf_core/lint.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir)
```

Initialise linting object

---

<a href="../../../../../../tools/nf_core/lint.py#L579"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_anaconda_package`

```python
check_anaconda_package(dep)
```

Call the anaconda API to find details about package

---

<a href="../../../../../../tools/nf_core/lint.py#L371"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_ci_config`

```python
check_ci_config()
```

Check that the Travis or Circle CI YAML config is valid

Makes sure that `nf-core lint` runs in travis tests Checks that tests run with the required nextflow version

---

<a href="../../../../../../tools/nf_core/lint.py#L623"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_conda_dockerfile`

```python
check_conda_dockerfile()
```

Check that the Docker build file looks right, if working with conda

Make sure that a name is given and is consistent with the pipeline name Check that depedency versions are pinned Warn if dependency versions are not the latest available

---

<a href="../../../../../../tools/nf_core/lint.py#L516"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_conda_env_yaml`

```python
check_conda_env_yaml()
```

Check that the conda environment file is valid

Make sure that a name is given and is consistent with the pipeline name Check that depedency versions are pinned Warn if dependency versions are not the latest available

---

<a href="../../../../../../tools/nf_core/lint.py#L647"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_conda_singularityfile`

```python
check_conda_singularityfile()
```

Check that the Singularity build file looks right, if working with conda

Make sure that a name is given and is consistent with the pipeline name Check that depedency versions are pinned Warn if dependency versions are not the latest available

---

<a href="../../../../../../tools/nf_core/lint.py#L194"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_docker`

```python
check_docker()
```

Check that Dockerfile contains the string 'FROM '

---

<a href="../../../../../../tools/nf_core/lint.py#L134"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_files_exist`

```python
check_files_exist()
```

Check a given pipeline directory for required files.

Throws an AssertionError if neither nextflow.config or main.nf found Gives either test failures or warnings for set of other filenames

---

<a href="../../../../../../tools/nf_core/lint.py#L223"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_licence`

```python
check_licence()
```

Check licence file is MIT

Ensures that Licence file is long enough (4 or more lines) Checks that licence contains the string 'without restriction' Checks that licence doesn't have any placeholder variables

---

<a href="../../../../../../tools/nf_core/lint.py#L267"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_nextflow_config`

```python
check_nextflow_config()
```

Check a given pipeline for required config variables.

Uses `nextflow config -flat` to parse pipeline nextflow.config and print all config variables. NB: Does NOT parse contents of main.nf / nextflow script

---

<a href="../../../../../../tools/nf_core/lint.py#L603"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_pip_package`

```python
check_pip_package(dep)
```

Call the PyPI API to find details about package

---

<a href="../../../../../../tools/nf_core/lint.py#L437"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_readme`

```python
check_readme()
```

Check the repository README file for errors

Currently just checks the badges at the top of the README

---

<a href="../../../../../../tools/nf_core/lint.py#L208"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_singularity`

```python
check_singularity()
```

Check that Singularity file contains the string 'FROM '

---

<a href="../../../../../../tools/nf_core/lint.py#L469"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_version_consistency`

```python
check_version_consistency()
```

Check container tags versions

Runs on process.container, params.container and $TRAVIS_TAG (each only if set) Check that the container has a tag Check that the version numbers are numeric Check that the version numbers are the same as one-another

---

<a href="../../../../../../tools/nf_core/lint.py#L84"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_pipeline`

```python
lint_pipeline(release=False)
```

Main linting function.

Takes the pipeline directory as the primary input and iterates through the different linting checks in order. Collects any warnings or errors and returns summary at completion. Raises an exception if there is a critical error that makes the rest of the tests pointless (eg. no pipeline script). Results from this function are printed by the main script.

**Args:**

- <b>`pipeline_dir`</b> (str): The path to the pipeline directory

**Returns:**

- <b>`dict`</b>: Summary of test result messages structured as follows: {
- <b>`'pass'`</b>: [ ( test-id (int), message (string) ), ( test-id (int), message (string) ) ],
- <b>`'warn'`</b>: [(id, msg)],
- <b>`'fail'`</b>: [(id, msg)], }

**Raises:**
If a critical problem is found, an AssertionError is raised.

---

<a href="../../../../../../tools/nf_core/lint.py#L675"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_results`

```python
print_results()
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
