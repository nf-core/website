<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint`

Code for linting modules in the nf-core/modules repository and in nf-core pipelines

Command: nf-core modules lint

## **Global Variables**

- **main_nf**
- **functions_nf**
- **meta_yml**
- **module_changes**
- **module_tests**
- **module_todos**
- **module_version**

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLintException`

Exception raised when there was an error with module linting

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L45"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `LintResult`

An object to hold the results of a lint test

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(mod, lint_test, message, file_path)
```

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L56"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLint`

An object for linting modules either in a clone of the 'nf-core/modules' repository or in any nf-core pipeline directory

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L71"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(dir)
```

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L192"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `filter_tests_by_key`

```python
filter_tests_by_key(key)
```

Filters the tests by the supplied key

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L207"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_installed_modules`

```python
get_installed_modules()
```

Makes lists of the local and and nf-core modules installed in this directory.

**Returns:**
local_modules, nfcore_modules ([NfCoreModule], [NfCoreModule]):

- <b>`A tuple of two lists`</b>: One for local modules and one for nf-core modules. In case the module contains several subtools, one path to each tool directory is returned.

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L98"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint`

```python
lint(
    module=None,
    key=(),
    all_modules=False,
    print_results=True,
    show_passed=False,
    local=False
)
```

Lint all or one specific module

First gets a list of all local modules (in modules/local/process) and all modules installed from nf-core (in modules/nf-core/modules)

For all nf-core modules, the correct file structure is assured and important file content is verified. If directory subject to linting is a clone of 'nf-core/modules', the files necessary for testing the modules are also inspected.

For all local modules, the '.nf' file is checked for some important flags, and warnings are issued if untypical content is found.

:param module: A specific module to lint :param print_results: Whether to print the linting results :param show_passed: Whether passed tests should be shown as well

:returns: A ModuleLint object containing information of the passed, warned and failed tests

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L299"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_module`

```python
lint_module(mod, local=False)
```

Perform linting on one module

If the module is a local module we only check the `main.nf` file, and issue warnings instead of failures.

If the module is a nf-core module we check for existence of the files

- main.nf
- meta.yml
- functions.nf And verify that their content conform to the nf-core standards.

If the linting is run for modules in the central nf-core/modules repo (repo_type==modules), files that are relevant for module testing are also examined

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L274"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_modules`

```python
lint_modules(modules, local=False)
```

Lint a list of modules

**Args:**

- <b>`modules`</b> ([NFCoreModule]): A list of module objects
- <b>`local`</b> (boolean): Whether the list consist of local or nf-core modules

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L426"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_summary`

```python
print_summary()
```

---

<a href="../../../../../../tools/nf_core/modules/lint/__init__.py#L181"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `set_up_pipeline_files`

```python
set_up_pipeline_files()
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
