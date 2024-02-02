# nf_core.bump_version

Bumps the version number in all appropriate files for
a nf-core pipeline.

### nf_core.bump_version.bump_nextflow_version(pipeline_obj: [Pipeline](utils.md#nf_core.utils.Pipeline), new_version: str)

Bumps the required Nextflow version number of a pipeline.

- **Parameters:**
  - **pipeline_obj** ([_nf_core.utils.Pipeline_](utils.md#nf_core.utils.Pipeline)) – A Pipeline object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the required Nextflow version.

### nf_core.bump_version.bump_pipeline_version(pipeline_obj: [Pipeline](utils.md#nf_core.utils.Pipeline), new_version: str)

Bumps a pipeline version number.

- **Parameters:**
  - **pipeline_obj** ([_nf_core.utils.Pipeline_](utils.md#nf_core.utils.Pipeline)) – A Pipeline object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the pipeline. Semantic versioning only.

### nf_core.bump_version.update_file_version(filename: str | Path, pipeline_obj: [Pipeline](utils.md#nf_core.utils.Pipeline), patterns: List[Tuple[str, str]])

Updates the version number in a requested file.

- **Parameters:**
  - **filename** (_str_) – File to scan.
  - **pipeline_obj** ([_nf_core.lint.PipelineLint_](lint.md#nf_core.lint.PipelineLint)) – A PipelineLint object that holds information
    about the pipeline contents and build files.
  - **pattern** (_str_) – Regex pattern to apply.
- **Raises:**
  **ValueError\*\***,\*\* **if the version number cannot be found.** –
