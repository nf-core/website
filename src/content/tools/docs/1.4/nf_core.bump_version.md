<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/bump_version.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.bump_version`

Release code for the nf-core python package.

Bumps the version number in all appropriate files for a nf-core pipeline

---

<a href="../../../../../../tools/nf_core/bump_version.py#L13"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `bump_pipeline_version`

```python
bump_pipeline_version(lint_obj, new_version)
```

Function to bump a pipeline version number. Called by the main script

---

<a href="../../../../../../tools/nf_core/bump_version.py#L67"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `bump_nextflow_version`

```python
bump_nextflow_version(lint_obj, new_version)
```

Function to bump the required nextflow version number.

---

<a href="../../../../../../tools/nf_core/bump_version.py#L94"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `update_file_version`

```python
update_file_version(filename, lint_obj, pattern, newstr, allow_multiple=False)
```

Update version number in the requested file

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
