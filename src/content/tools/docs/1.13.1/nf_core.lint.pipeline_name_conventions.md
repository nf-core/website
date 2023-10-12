<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/pipeline_name_conventions.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.pipeline_name_conventions`

---

<a href="../../../../../../tools/nf_core/lint/pipeline_name_conventions.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `pipeline_name_conventions`

```python
pipeline_name_conventions()
```

Checks that the pipeline name adheres to nf-core conventions.

In order to ensure consistent naming, pipeline names should contain only lower case, alphanumeric characters. Otherwise a warning is displayed.

.. warning:
`` DockerHub is very picky about image names and doesn't even allow hyphens (we are`nfcore``). This is a large part of why we set this rule.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
