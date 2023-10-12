<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/system_exit.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.system_exit`

---

<a href="../../../../../../tools/nf_core/lint/system_exit.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `system_exit`

```python
system_exit()
```

Check for System.exit calls in groovy/nextflow code

Calls to System.exit(1) should be replaced by throwing errors

This lint test looks for all calls to `System.exit` in any file with the `.nf` or `.groovy` extension

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
