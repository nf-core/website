<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/bump_version.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.bump_version`

Bumps the version number in all appropriate files for a nf-core pipeline.

---

<a href="../../../../../../tools/nf_core/bump_version.py#L18"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `bump_pipeline_version`

```python
bump_pipeline_version(pipeline_obj, new_version)
```

Bumps a pipeline version number.

**Args:**

- <b>`pipeline_obj`</b> (nf_core.utils.Pipeline): A `Pipeline` object that holds information about the pipeline contents and build files.
- <b>`new_version`</b> (str): The new version tag for the pipeline. Semantic versioning only.

---

<a href="../../../../../../tools/nf_core/bump_version.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `bump_nextflow_version`

```python
bump_nextflow_version(pipeline_obj, new_version)
```

Bumps the required Nextflow version number of a pipeline.

**Args:**

- <b>`pipeline_obj`</b> (nf_core.utils.Pipeline): A `Pipeline` object that holds information about the pipeline contents and build files.
- <b>`new_version`</b> (str): The new version tag for the required Nextflow version.

---

<a href="../../../../../../tools/nf_core/bump_version.py#L167"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `update_file_version`

```python
update_file_version(filename, pipeline_obj, patterns)
```

Updates the version number in a requested file.

**Args:**

- <b>`filename`</b> (str): File to scan.
- <b>`pipeline_obj`</b> (nf_core.lint.PipelineLint): A PipelineLint object that holds information about the pipeline contents and build files.
- <b>`pattern`</b> (str): Regex pattern to apply.
- <b>`newstr`</b> (str): The replaced string.

**Raises:**
ValueError, if the version number cannot be found.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
