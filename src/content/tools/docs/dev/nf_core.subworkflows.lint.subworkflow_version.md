<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_version.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint.subworkflow_version`

Verify that a subworkflow has a correct entry in the modules.json file

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_version.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `subworkflow_version`

```python
subworkflow_version(subworkflow_lint_object, subworkflow)
```

Verifies that the subworkflow has a version specified in the `modules.json` file

It checks whether the subworkflow has an entry in the `modules.json` file containing a commit SHA. If that is true, it verifies that there are no newer version of the subworkflow available.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
