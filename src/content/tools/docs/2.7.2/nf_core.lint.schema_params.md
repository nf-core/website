<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/schema_params.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.schema_params`

---

<a href="../../../../../../tools/nf_core/lint/schema_params.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `schema_params`

```python
schema_params()
```

Check that the schema describes all flat params in the pipeline.

The `nextflow_schema.json` pipeline schema should describe every flat parameter returned from the `nextflow config` command (params that are objects or more complex structures are ignored).

- Failure: If parameters are found in `nextflow_schema.json` that are not in `nextflow_schema.json` \* Warning: If parameters are found in `nextflow_schema.json` that are not in `nextflow_schema.json`

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
