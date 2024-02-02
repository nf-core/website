# nf_core.bump_version

Bumps the version number in all appropriate files for
a nf-core pipeline.

### nf_core.bump_version.bump_nextflow_version(lint_obj, new_version)

Bumps the required Nextflow version number of a pipeline.

- **Parameters:**
  - **lint_obj** ([_nf_core.lint.PipelineLint_](lint.md#nf_core.lint.PipelineLint)) – A PipelineLint object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the required Nextflow version.

### nf_core.bump_version.bump_pipeline_version(lint_obj, new_version)

Bumps a pipeline version number.

- **Parameters:**
  - **lint_obj** ([_nf_core.lint.PipelineLint_](lint.md#nf_core.lint.PipelineLint)) – A PipelineLint object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the pipeline. Semantic versioning only.

### nf_core.bump_version.update_file_version(filename, lint_obj, pattern, newstr, allow_multiple=False)

Updates the version number in a requested file.

- **Parameters:**
  - **filename** (_str_) – File to scan.
  - **lint_obj** ([_nf_core.lint.PipelineLint_](lint.md#nf_core.lint.PipelineLint)) – A PipelineLint object that holds information
    about the pipeline contents and build files.
  - **pattern** (_str_) – Regex pattern to apply.
  - **newstr** (_str_) – The replaced string.
  - **allow_multiple** (_bool_) – Replace all pattern hits, not only the first. Defaults to False.
- **Raises:**
  **SyntaxError\*\***,\*\* **if the version number cannot be found.** –
