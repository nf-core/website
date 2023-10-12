<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_changes.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint.subworkflow_changes`

Check whether the content of a subworkflow has changed compared to the original repository

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_changes.py#L9"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `subworkflow_changes`

```python
subworkflow_changes(subworkflow_lint_object, subworkflow)
```

Checks whether installed nf-core subworkflow have changed compared to the original repository

Downloads the `main.nf` and `meta.yml` files for every subworkflow and compares them to the local copies

If the subworkflow has a commit SHA entry in the `modules.json`, the file content is compared against the files in the remote at this SHA.

Only runs when linting a pipeline, not the modules repository

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
