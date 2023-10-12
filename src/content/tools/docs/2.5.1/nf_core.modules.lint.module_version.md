<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/module_version.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.module_version`

Verify that a module has a correct entry in the modules.json file

---

<a href="../../../../../../tools/nf_core/modules/lint/module_version.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `module_version`

```python
module_version(module_lint_object, module)
```

Verifies that the module has a version specified in the `modules.json` file

It checks whether the module has an entry in the `modules.json` file containing a commit SHA. If that is true, it verifies that there are no newer version of the module available.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
