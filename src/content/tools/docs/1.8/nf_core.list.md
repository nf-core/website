<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/list.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.list`

Lists available nf-core pipelines and versions.

---

<a href="../../../../../../tools/nf_core/list.py#L27"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `list_workflows`

```python
list_workflows(filter_by=None, sort_by='release', as_json=False)
```

Prints out a list of all nf-core workflows.

**Args:**

- <b>`filter_by`</b> (list): A list of strings that can be used for filtering.
- <b>`sort_by`</b> (str): workflows can be sorted by keywords. Keyword must be one of `release` (default), `name`, `stars`.
- <b>`as_json`</b> (boolean): Set to true, if the lists should be printed in JSON.

---

<a href="../../../../../../tools/nf_core/list.py#L325"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `pretty_date`

```python
pretty_date(time)
```

Transforms a datetime object or a int() Epoch timestamp into a pretty string like 'an hour ago', 'Yesterday', '3 months ago', 'just now', etc

Based on https://stackoverflow.com/a/1551394/713980 Adapted by sven1103

---

<a href="../../../../../../tools/nf_core/list.py#L46"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Workflows`

Workflow container class.

Is used to collect local and remote nf-core pipelines. Pipelines can be sorted, filtered and compared.

**Args:**

- <b>`filter_by`</b> (list): A list of strings that can be used for filtering.
- <b>`sort_by`</b> (str): workflows can be sorted by keywords. Keyword must be one of `release` (default), `name`, `stars`.

<a href="../../../../../../tools/nf_core/list.py#L57"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(filter_by=None, sort_by='release')
```

---

<a href="../../../../../../tools/nf_core/list.py#L116"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compare_remote_local`

```python
compare_remote_local()
```

Matches local to remote workflows.

If a matching remote workflow is found, the local workflow's Git commit hash is compared with the latest one from remote.

A boolean flag in :attr:`RemoteWorkflow.local_is_latest` is set to True, if the local workflow is the latest.

---

<a href="../../../../../../tools/nf_core/list.py#L135"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `filtered_workflows`

```python
filtered_workflows()
```

Filters remote workflows for keywords.

**Returns:**

- <b>`list`</b>: Filtered remote workflows.

---

<a href="../../../../../../tools/nf_core/list.py#L78"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_nf_workflows`

```python
get_local_nf_workflows()
```

Retrieves local Nextflow workflows.

Local workflows are stored in :attr:`self.local_workflows` list.

---

<a href="../../../../../../tools/nf_core/list.py#L64"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_workflows`

```python
get_remote_workflows()
```

Retrieves remote workflows from `nf-co.re <http://nf-co.re>`\_.

Remote workflows are stored in :attr:`self.remote_workflows` list.

---

<a href="../../../../../../tools/nf_core/list.py#L214"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_json`

```python
print_json()
```

Dump JSON of all parsed information

---

<a href="../../../../../../tools/nf_core/list.py#L158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_summary`

```python
print_summary()
```

Prints a summary of all pipelines.

---

<a href="../../../../../../tools/nf_core/list.py#L222"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `RemoteWorkflow`

A information container for a remote workflow.

**Args:**

- <b>`data`</b> (dict): workflow information as they are retrieved from the Github repository REST API request
- <b>`(https`</b>: //developer.github.com/v3/repos/#get).

<a href="../../../../../../tools/nf_core/list.py#L230"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(data)
```

---

<a href="../../../../../../tools/nf_core/list.py#L256"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `LocalWorkflow`

Class to handle local workflows pulled by nextflow

<a href="../../../../../../tools/nf_core/list.py#L259"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(name)
```

Initialise the LocalWorkflow object

---

<a href="../../../../../../tools/nf_core/list.py#L271"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_nf_workflow_details`

```python
get_local_nf_workflow_details()
```

Get full details about a local cached workflow

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
