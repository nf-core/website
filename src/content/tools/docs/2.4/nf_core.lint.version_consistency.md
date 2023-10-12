<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/version_consistency.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.version_consistency`

---

<a href="../../../../../../tools/nf_core/lint/version_consistency.py#L6"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `version_consistency`

```python
version_consistency()
```

Pipeline and container version number consistency.

.. note:: This test only runs when the `--release` flag is set for `nf-core lint`, or `$GITHUB_REF` is equal to `master`.

This lint fetches the pipeline version number from three possible locations:

- The pipeline config, `manifest.version` \* The docker container in the pipeline config, `process.container`

- Some pipelines may not have this set on a pipeline level. If it is not found, it is ignored.

- `$GITHUB_REF`, if it looks like a release tag (`refs/tags/<something>`)

The test then checks that:

- The container name has a tag specified (eg. `nfcore/pipeline:version`) _ The pipeline version number is numeric (contains only numbers and dots) _ That the version numbers all match one another

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
