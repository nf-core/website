<!-- markdownlint-disable -->

<a href="../../nf_core/lint/actions_schema_validation.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.actions_schema_validation`





---

<a href="../../nf_core/lint/actions_schema_validation.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `actions_schema_validation`

```python
actions_schema_validation()
```

Checks that the GitHub Action workflow yml/yaml files adhere to the correct schema 

nf-core pipelines use GitHub actions workflows to run CI tests, check formatting and also linting, among others. These workflows are defined by ``yml`` scripts in ``.github/workflows/``. This lint test verifies that these scripts are valid by comparing them against the `JSON schema for GitHub workflows <https://json.schemastore.org/github-workflow>`_. 

To pass this test, make sure that all your workflows contain the required properties ``on`` and ``jobs`` and that all other properties are of the correct type, as specified in the schema (link above). 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
