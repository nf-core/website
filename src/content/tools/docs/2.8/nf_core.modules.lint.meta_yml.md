<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/meta_yml.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.meta_yml`

---

<a href="../../../../../../tools/nf_core/modules/lint/meta_yml.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `meta_yml`

```python
meta_yml(module_lint_object, module)
```

Lint a `meta.yml` file

The lint test checks that the module has a `meta.yml` file and that it follows the JSON schema defined in the `modules/yaml-schema.json` file in the nf-core/modules repository.

In addition it checks that the module name and module input is consistent between the `meta.yml` and the `main.nf`.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
