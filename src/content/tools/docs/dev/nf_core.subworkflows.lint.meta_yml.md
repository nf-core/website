<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/meta_yml.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint.meta_yml`

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/meta_yml.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `meta_yml`

```python
meta_yml(subworkflow_lint_object, subworkflow)
```

Lint a `meta.yml` file

The lint test checks that the subworkflow has a `meta.yml` file and that it follows the JSON schema defined in the `subworkflows/yaml-schema.json` file in the nf-core/modules repository.

In addition it checks that the subworkflow name and subworkflow input is consistent between the `meta.yml` and the `main.nf`.

Checks that all input and output channels are specified in `meta.yml`. Checks that all included components in `main.nf` are specified in `meta.yml`.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
