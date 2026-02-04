# nf_core.pipelines.lint

:::tip
See the [Lint Tests]() docs for information about specific linting functions.
:::

<a id="module-nf_core.pipelines.lint"></a>

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to
the nf-core community guidelines.

### `nf_core.pipelines.lint.run_linting(pipeline_dir, release_mode: bool = False, fix=(), key=(), show_passed: bool = False, fail_ignored: bool = False, fail_warned: bool = False, sort_by: str = 'test', md_fn=None, json_fn=None, hide_progress: bool = False, plain_text: bool = False) → tuple[PipelineLint, ComponentLint | None, ComponentLint | None]{:python}`

Runs all nf-core linting checks on a given Nextflow pipeline project
in either release mode or normal mode (default). Returns an object
of type `PipelineLint` after finished.

- **Parameters:**
  - **pipeline_dir** (_str_) – The path to the Nextflow pipeline root directory
  - **release_mode** (_bool_) – Set this to True, if the linting should be run in the release mode.
    See `PipelineLint` for more information.
  - **plain_text** (_bool_) – Print output in plain text without rich formatting
- **Returns:**
  An object of type `PipelineLint` that contains all the linting results.
  An object of type `ComponentLint` that contains all the linting results for the modules.
  An object of type `ComponentLint` that contains all the linting results for the subworkflows.
