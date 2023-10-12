<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/list.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.list`

List available nf-core pipelines and versions

---

<a href="../../../../../../tools/nf_core/list.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `list_workflows`

```python
list_workflows(sort='release', json=False, keywords=[])
```

Main function to list all nf-core workflows

---

<a href="../../../../../../tools/nf_core/list.py#L268"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `pretty_date`

```python
pretty_date(time)
```

Get a datetime object or a int() Epoch timestamp and return a pretty string like 'an hour ago', 'Yesterday', '3 months ago', 'just now', etc

Based on https://stackoverflow.com/a/1551394/713980 Adapted by sven1103

---

<a href="../../../../../../tools/nf_core/list.py#L35"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Workflows`

Class to hold all workflows

<a href="../../../../../../tools/nf_core/list.py#L38"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(sort='release', keywords=[])
```

Initialise the class with empty placeholder vars

---

<a href="../../../../../../tools/nf_core/list.py#L94"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `compare_remote_local`

```python
compare_remote_local()
```

Match local to remote workflows.

---

<a href="../../../../../../tools/nf_core/list.py#L106"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `filtered_workflows`

```python
filtered_workflows()
```

Filter remote workflows if keywords supplied

---

<a href="../../../../../../tools/nf_core/list.py#L58"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_nf_workflows`

```python
get_local_nf_workflows()
```

Get local nextflow workflows

---

<a href="../../../../../../tools/nf_core/list.py#L46"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_remote_workflows`

```python
get_remote_workflows()
```

Get remote nf-core workflows

---

<a href="../../../../../../tools/nf_core/list.py#L170"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_json`

```python
print_json()
```

Dump JSON of all parsed information

---

<a href="../../../../../../tools/nf_core/list.py#L125"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_summary`

```python
print_summary()
```

Print summary of all pipelines

---

<a href="../../../../../../tools/nf_core/list.py#L178"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `RemoteWorkflow`

Class to hold a single workflow

<a href="../../../../../../tools/nf_core/list.py#L181"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(data)
```

Initialise a workflow object from the GitHub API object

---

<a href="../../../../../../tools/nf_core/list.py#L209"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `LocalWorkflow`

Class to handle local workflows pulled by nextflow

<a href="../../../../../../tools/nf_core/list.py#L212"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(name)
```

Initialise the LocalWorkflow object

---

<a href="../../../../../../tools/nf_core/list.py#L224"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_local_nf_workflow_details`

```python
get_local_nf_workflow_details()
```

Get full details about a local cached workflow

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
