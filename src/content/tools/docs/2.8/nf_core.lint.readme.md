<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/readme.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.readme`

---

<a href="../../../../../../tools/nf_core/lint/readme.py#L5"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `readme`

```python
readme()
```

Repository `README.md` tests

The `README.md` files for a project are very important and must meet some requirements:

- Nextflow badge

- If no Nextflow badge is found, a warning is given _ If a badge is found but the version doesn't match the minimum version in the config file, the test fails _ Example badge code:

.. code-block:: md

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)

.. note:: This badge are a markdown image `![alt-text](<image URL>)` _inside_ a markdown link `[markdown image](<link URL>)`, so a bit fiddly to write.

- Zenodo release

- If pipeline is released but still contains a 'zenodo.XXXXXXX' tag, the test fails

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
