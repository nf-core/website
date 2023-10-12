<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/licences.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.licences`

Lists software licences for a given workflow.

---

<a href="../../../../../../tools/nf_core/licences.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `WorkflowLicences`

A nf-core workflow licenses collection.

Tries to retrieve the license information from all dependencies of a given nf-core pipeline.

A condensed overview with license per dependency can be printed out.

**Args:**

- <b>`pipeline`</b> (str): An existing nf-core pipeline name, like `nf-core/hlatyping` or short `hlatyping`.

<a href="../../../../../../tools/nf_core/licences.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline)
```

---

<a href="../../../../../../tools/nf_core/licences.py#L71"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `fetch_conda_licences`

```python
fetch_conda_licences()
```

Fetch package licences from Anaconda and PyPi.

---

<a href="../../../../../../tools/nf_core/licences.py#L50"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_environment_file`

```python
get_environment_file()
```

Get the conda environment file for the pipeline

---

<a href="../../../../../../tools/nf_core/licences.py#L92"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_licences`

```python
print_licences()
```

Prints the fetched license information.

**Args:**

- <b>`as_json`</b> (boolean): Prints the information in JSON. Defaults to False.

---

<a href="../../../../../../tools/nf_core/licences.py#L42"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `run_licences`

```python
run_licences()
```

Run the nf-core licences action

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
