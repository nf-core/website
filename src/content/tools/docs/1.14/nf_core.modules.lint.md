<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint`

Code for linting modules in the nf-core/modules repository and in nf-core pipelines

Command: nf-core modules lint

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLintException`

Exception raised when there was an error with module linting

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `LintResult`

An object to hold the results of a lint test

<a href="../../../../../../tools/nf_core/modules/lint.py#L42"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(mod, lint_test, message, file_path)
```

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L50"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLint`

An object for linting modules either in a clone of the 'nf-core/modules' repository or in any nf-core pipeline directory

<a href="../../../../../../tools/nf_core/modules/lint.py#L56"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(dir)
```

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L373"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_module_changes`

```python
check_module_changes(nfcore_modules)
```

Checks whether installed nf-core modules have changed compared to the original repository Downloads the 'main.nf', 'functions.nf' and 'meta.yml' files for every module and compare them to the local copies

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L214"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_installed_modules`

```python
get_installed_modules()
```

Make a list of all modules installed in this repository

Returns a tuple of two lists, one for local modules and one for nf-core modules. The local modules are represented as direct filepaths to the module '.nf' file. Nf-core module are returned as file paths to the module directories. In case the module contains several tools, one path to each tool directory is returned.

returns (local_modules, nfcore_modules)

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L195"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_repo_type`

```python
get_repo_type()
```

Determine whether this is a pipeline repository or a clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L64"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint`

```python
lint(
    module=None,
    all_modules=False,
    print_results=True,
    show_passed=False,
    local=False
)
```

Lint all or one specific module

First gets a list of all local modules (in modules/local/process) and all modules installed from nf-core (in modules/nf-core/software)

For all nf-core modules, the correct file structure is assured and important file content is verified. If directory subject to linting is a clone of 'nf-core/modules', the files necessary for testing the modules are also inspected.

For all local modules, the '.nf' file is checked for some important flags, and warnings are issued if untypical content is found.

:param module: A specific module to lint :param print_results: Whether to print the linting results :param show_passed: Whether passed tests should be shown as well

:returns: dict of {passed, warned, failed}

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L137"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_local_modules`

```python
lint_local_modules(local_modules)
```

Lint a local module Only issues warnings instead of failures

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L164"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_nfcore_modules`

```python
lint_nfcore_modules(nfcore_modules)
```

Lint nf-core modules For each nf-core module, checks for existence of the files

- main.nf
- meta.yml
- functions.nf And verifies that their content.

If the linting is run for modules in the central nf-core/modules repo (repo_type==modules), files that are relevant for module testing are also examined

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L450"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `NFCoreModule`

A class to hold the information a bout a nf-core module Includes functionality for linting

<a href="../../../../../../tools/nf_core/modules/lint.py#L456"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(module_dir, repo_type, base_dir, nf_core_module=True)
```

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L730"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_process_section`

```python
check_process_section(lines)
```

Lint the section of a module between the process definition and the 'input:' definition Specifically checks for correct software versions and containers

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L708"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_script_section`

```python
check_script_section(lines)
```

Lint the script section Checks whether 'def sotware' and 'def prefix' are defined

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L478"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint`

```python
lint()
```

Perform linting on this module

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L811"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_functions_nf`

```python
lint_functions_nf()
```

Lint a functions.nf file Verifies that the file exists and contains all necessary functions

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L603"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_main_nf`

```python
lint_main_nf()
```

Lint a single main.nf module file Can also be used to lint local module files, in which case failures should be interpreted as warnings

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L554"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_meta_yml`

```python
lint_meta_yml()
```

Lint a meta yml file

---

<a href="../../../../../../tools/nf_core/modules/lint.py#L503"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_module_tests`

```python
lint_module_tests()
```

Lint module tests

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
