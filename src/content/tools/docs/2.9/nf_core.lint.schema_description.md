<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/schema_description.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.schema_description`

---

<a href="../../../../../../tools/nf_core/lint/schema_description.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `schema_description`

```python
schema_description()
```

Check that every parameter in the schema has a description

The `nextflow_schema.json` pipeline schema should describe every flat parameter Furthermore warns about parameters outside of groups

- Warning: Parameters in `nextflow_schema.json` without a description \* Warning: Parameters in `nextflow_schema.json` that are defined outside of a group

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
