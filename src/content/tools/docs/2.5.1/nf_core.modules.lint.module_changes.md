<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/module_changes.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.module_changes`

Check whether the content of a module has changed compared to the original repository

---

<a href="../../../../../../tools/nf_core/modules/lint/module_changes.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `module_changes`

```python
module_changes(module_lint_object, module)
```

Checks whether installed nf-core modules have changed compared to the original repository

Downloads the `main.nf` and `meta.yml` files for every module and compares them to the local copies

If the module has a commit SHA entry in the `modules.json`, the file content is compared against the files in the remote at this SHA.

Only runs when linting a pipeline, not the modules repository

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
