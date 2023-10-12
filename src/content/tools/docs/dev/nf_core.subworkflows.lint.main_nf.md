<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/main_nf.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint.main_nf`

Lint the main.nf file of a subworkflow

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/main_nf.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main_nf`

```python
main_nf(_, subworkflow)
```

Lint a `main.nf` subworkflow file

Can also be used to lint local subworkflow files, in which case failures will be reported as warnings.

The test checks for the following:

- A subworkflow SHOULD import at least two modules _ All included modules or subworkflows are used and their names are used for `versions.yml` _ The workflow name is all capital letters \* The subworkflow emits a software version

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/main_nf.py#L102"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_main_section`

```python
check_main_section(lines, included_components)
```

Lint the main section of a subworkflow Checks whether all included components are used and their names are used for `versions.yml`.

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/main_nf.py#L155"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_subworkflow_section`

```python
check_subworkflow_section(lines)
```

Lint the section of a subworkflow before the workflow definition Specifically checks if the subworkflow includes at least two modules or subworkflows

**Args:**

- <b>`lines`</b> (List[str]): Content of subworkflow.

**Returns:**

- <b>`List`</b>: List of included component names. If subworkflow doesn't contain any lines, return None.

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/main_nf.py#L194"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_workflow_section`

```python
check_workflow_section(lines)
```

Lint the workflow definition of a subworkflow before Specifically checks that the name is all capital letters

**Args:**

- <b>`lines`</b> (List[str]): Content of workflow definition.

**Returns:**
None

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
