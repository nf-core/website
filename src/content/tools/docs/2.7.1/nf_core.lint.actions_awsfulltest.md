<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/actions_awsfulltest.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.actions_awsfulltest`

---

<a href="../../../../../../tools/nf_core/lint/actions_awsfulltest.py#L6"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `actions_awsfulltest`

```python
actions_awsfulltest()
```

Checks the GitHub Actions awsfulltest is valid.

In addition to small test datasets run on GitHub Actions, we provide the possibility of testing the pipeline on full size datasets on AWS. This should ensure that the pipeline runs as expected on AWS and provide a resource estimation.

The GitHub Actions workflow is called `awsfulltest.yml`, and it can be found in the `.github/workflows/` directory.

.. warning:: This workflow incurs AWS costs, therefore it should only be triggered for pipeline releases: `release` (after the pipeline release) and `workflow_dispatch`.

.. note:: You can manually trigger the AWS tests by going to the `Actions` tab on the pipeline GitHub repository and selecting the `nf-core AWS full size tests` workflow on the left.

.. tip:: For tests on full data prior to release, `Nextflow Tower <https://tower.nf>`\_ launch feature can be employed.

The `.github/workflows/awsfulltest.yml` file is tested for the following:

- Must be turned on `workflow_dispatch`. _ Must be turned on for `release` with `types: [published]` _ Should run the profile `test_full` that should be edited to provide the links to full-size datasets. If it runs the profile `test`, a warning is given.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
