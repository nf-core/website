<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.lint`

Code for linting subworkflows in the nf-core/subworkflows repository and in nf-core pipelines

Command: nf-core subworkflows lint

## **Global Variables**

- **main_nf**
- **meta_yml**
- **subworkflow_changes**
- **subworkflow_tests**
- **subworkflow_todos**
- **subworkflow_version**

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L25"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SubworkflowLint`

An object for linting subworkflows either in a clone of the 'nf-core/modules' repository or in any nf-core pipeline directory

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    dir,
    fail_warned=False,
    remote_url=None,
    branch=None,
    no_pull=False,
    registry=None,
    hide_progress=False
)
```

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L60"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint`

```python
lint(
    subworkflow=None,
    registry='quay.io',
    key=(),
    all_subworkflows=False,
    print_results=True,
    show_passed=False,
    sort_by='test',
    local=False
)
```

Lint all or one specific subworkflow

First gets a list of all local subworkflows (in subworkflows/local/process) and all subworkflows installed from nf-core (in subworkflows/nf-core)

For all nf-core subworkflows, the correct file structure is assured and important file content is verified. If directory subject to linting is a clone of 'nf-core/modules', the files necessary for testing the subworkflows are also inspected.

For all local subworkflows, the '.nf' file is checked for some important flags, and warnings are issued if untypical content is found.

:param subworkflow: A specific subworkflow to lint :param print_results: Whether to print the linting results :param show_passed: Whether passed tests should be shown as well :param hide_progress: Don't show progress bars

:returns: A SubworkflowLint object containing information of the passed, warned and failed tests

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L183"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_subworkflow`

```python
lint_subworkflow(swf, progress_bar, registry, local=False)
```

Perform linting on one subworkflow

If the subworkflow is a local subworkflow we only check the `main.nf` file, and issue warnings instead of failures.

If the subworkflow is a nf-core subworkflow we check for existence of the files

- main.nf
- meta.yml And verify that their content conform to the nf-core standards.

If the linting is run for subworkflows in the central nf-core/modules repo (repo_type==modules), files that are relevant for subworkflow testing are also examined

---

<a href="../../../../../../tools/nf_core/subworkflows/lint/__init__.py#L154"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_subworkflows`

```python
lint_subworkflows(subworkflows, registry='quay.io', local=False)
```

Lint a list of subworkflows

**Args:**

- <b>`subworkflows`</b> ([NFCoreComponent]): A list of subworkflow objects
- <b>`registry`</b> (str): The container registry to use. Should be quay.io in most situations.
- <b>`local`</b> (boolean): Whether the list consist of local or nf-core subworkflows

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
