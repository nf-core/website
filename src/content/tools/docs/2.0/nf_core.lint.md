<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/__init__.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint`

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to the nf-core community guidelines.

## **Global Variables**

- **pipeline_todos**
- **actions_awsfulltest**
- **actions_awstest**
- **actions_ci**
- **actions_schema_validation**
- **files_exist**
- **files_unchanged**
- **merge_markers**
- **modules_json**
- **nextflow_config**
- **pipeline_name_conventions**
- **readme**
- **schema_lint**
- **schema_params**
- **schema_description**
- **template_strings**
- **version_consistency**

---

<a href="../../../../../../tools/nf_core/lint/__init__.py#L29"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `run_linting`

```python
run_linting(
    pipeline_dir,
    release_mode=False,
    fix=(),
    key=(),
    show_passed=False,
    fail_ignored=False,
    md_fn=None,
    json_fn=None
)
```

Runs all nf-core linting checks on a given Nextflow pipeline project in either `release` mode or `normal` mode (default). Returns an object of type :class:`PipelineLint` after finished.

**Args:**

- <b>`pipeline_dir`</b> (str): The path to the Nextflow pipeline root directory
- <b>`release_mode`</b> (bool): Set this to `True`, if the linting should be run in the `release` mode.
- <b>`See `</b>: class:`PipelineLint` for more information.

**Returns:**

- <b>`An object of type `</b>: class:`PipelineLint` that contains all the linting results.

---

<a href="../../../../../../tools/nf_core/lint/__init__.py#L102"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineLint`

Object to hold linting information and results.

Inherits :class:`nf_core.utils.Pipeline` class.

Use the :func:`PipelineLint._lint_pipeline` function to run lint tests.

**Args:**

- <b>`path`</b> (str): The path to the nf-core pipeline directory.

**Attributes:**

- <b>`failed`</b> (list): A list of tuples of the form: `(<test-name>, <reason>)`
- <b>`ignored`</b> (list): A list of tuples of the form: `(<test-name>, <reason>)`
- <b>`lint_config`</b> (dict): The parsed nf-core linting config for this pipeline
- <b>`passed`</b> (list): A list of tuples of the form: `(<test-name>, <reason>)`
- <b>`release_mode`</b> (bool): `True`, if you the to linting was run in release mode, `False` else.
- <b>`warned`</b> (list): A list of tuples of the form: `(<warned no>, <reason>)`

<a href="../../../../../../tools/nf_core/lint/__init__.py#L139"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(wf_path, release_mode=False, fix=(), key=(), fail_ignored=False)
```

Initialise linting object

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
