<!-- markdownlint-disable -->

<a href="../../nf_core/lint/actions_ci.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.actions_ci`





---

<a href="../../nf_core/lint/actions_ci.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `actions_ci`

```python
actions_ci()
```

Checks that the GitHub Actions pipeline CI (Continuous Integration) workflow is valid. 

The ``.github/workflows/ci.yml`` GitHub Actions workflow runs the pipeline on a minimal test dataset using ``-profile test`` to check that no breaking changes have been introduced. Final result files are not checked, just that the pipeline exists successfully. 

This lint test checks this GitHub Actions workflow file for the following: 

* Workflow must be triggered on the following events: 

 .. code-block:: yaml 

 on:  push:  branches: 
             - dev  pull_request:  release:  types: [published] 

* The minimum Nextflow version specified in the pipeline's ``nextflow.config`` matches that defined by ``nxf_ver`` in the test matrix: 

 .. code-block:: yaml  :emphasize-lines: 4 

 strategy:  matrix:  # Nextflow versions: check pipeline minimum and current latest  nxf_ver: ['19.10.0', ''] 

 .. note:: These ``matrix`` variables run the test workflow twice, varying the ``nxf_ver`` variable each time.  This is used in the ``nextflow run`` commands to test the pipeline with both the latest available version  of the pipeline (``''``) and the stated minimum required version. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
