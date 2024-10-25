# nf_core.lint

:::tip
See the [Lint Tests](../pipeline_lint_tests/index) docs for information about specific linting functions.
:::

<a id="module-nf_core.lint"></a>

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to
the nf-core community guidelines.

### `nf_core.lint.run_linting(pipeline_dir, release_mode: bool = False, fix=(), key=(), show_passed: bool = False, fail_ignored: bool = False, fail_warned: bool = False, sort_by: str = 'test', md_fn=None, json_fn=None, hide_progress: bool = False){:python}`

Runs all nf-core linting checks on a given Nextflow pipeline project
in either release mode or normal mode (default). Returns an object
of type [`PipelineLint`](#nf_core.lint.PipelineLint) after finished.

- **Parameters:**
  - **pipeline_dir** (_str_) – The path to the Nextflow pipeline root directory
  - **release_mode** (_bool_) – Set this to True, if the linting should be run in the release mode.
    See [`PipelineLint`](#nf_core.lint.PipelineLint) for more information.
- **Returns:**
  An object of type [`PipelineLint`](#nf_core.lint.PipelineLint) that contains all the linting results.
  An object of type `ComponentLint` that contains all the linting results for the modules.
  An object of type `ComponentLint` that contains all the linting results for the subworkflows.

### _`class{:python}`_`nf_core.lint.PipelineLint(wf_path, release_mode=False, fix=(), key=None, fail_ignored=False, fail_warned=False, hide_progress=False){:python}`

Bases: [`Pipeline`](utils#nf_core.utils.Pipeline)

Object to hold linting information and results.

Inherits [`nf_core.utils.Pipeline`](utils#nf_core.utils.Pipeline) class.

Use the [`PipelineLint._lint_pipeline()`](#nf_core.lint.PipelineLint._lint_pipeline) function to run lint tests.

- **Parameters:**
  **path** (_str_) – The path to the nf-core pipeline directory.

#### `failed{:python}`

A list of tuples of the form: `(<test-name>, <reason>)`

- **Type:**
  list

#### `ignored{:python}`

A list of tuples of the form: `(<test-name>, <reason>)`

- **Type:**
  list

#### `lint_config{:python}`

The parsed nf-core linting config for this pipeline

- **Type:**
  dict

#### `passed{:python}`

A list of tuples of the form: `(<test-name>, <reason>)`

- **Type:**
  list

#### `release_mode{:python}`

True, if you the to linting was run in release mode, False else.

- **Type:**
  bool

#### `warned{:python}`

A list of tuples of the form: `(<warned no>, <reason>)`

- **Type:**
  list

#### `_get_results_md(){:python}`

Create a markdown file suitable for posting in a GitHub comment.

- **Returns:**
  Formatting markdown content
- **Return type:**
  markdown (str)

#### `_lint_pipeline(){:python}`

Main linting function.

Takes the pipeline directory as the primary input and iterates through
the different linting checks in order. Collects any warnings or errors
into object attributes: `passed`, `ignored`, `warned` and `failed`.

#### `_print_results(show_passed){:python}`

Print linting results to the command line.

Uses the `rich` library to print a set of formatted tables to the command line
summarising the linting results.

#### `_save_json_results(json_fn){:python}`

Function to dump lint results to a JSON file for downstream use

- **Parameters:**
  **json_fn** (_str_) – File path to write JSON to.

#### `_wrap_quotes(files: List[str] | List[Path] | Path){:python}`

Helper function to take a list of filenames and format with markdown.

- **Parameters:**
  **files** (_list_) –

  List of filenames, eg:

  ```default
  ['foo', 'bar', 'baz']
  ```

- **Returns:**
  Formatted string of paths separated by word `or`, eg:

  ```default
  `foo` or bar` or `baz`
  ```

- **Return type:**
  markdown (str)
