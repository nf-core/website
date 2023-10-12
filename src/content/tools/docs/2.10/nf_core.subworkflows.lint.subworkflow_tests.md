<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_tests.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint.subworkflow_tests`

Lint the tests of a subworkflow in nf-core/modules

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/subworkflow_tests.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `subworkflow_tests`

```python
subworkflow_tests(_, subworkflow)
```

Lint the tests of a subworkflow in `nf-core/modules`

It verifies that the test directory exists and contains a `main.nf` and a `test.yml`, and that the subworkflow is present in the `pytest_modules.yml` file.

Additionally, hecks that all included components in test `main.nf` are specified in `test.yml`

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
