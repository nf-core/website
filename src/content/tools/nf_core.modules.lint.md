<!-- markdownlint-disable -->

<a href="../../nf_core/modules/lint/__init__.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint`
Code for linting modules in the nf-core/modules repository and in nf-core pipelines 

Command: nf-core modules lint 

**Global Variables**
---------------
- **main_nf**
- **meta_yml**
- **module_changes**
- **module_deprecations**
- **module_patch**
- **module_tests**
- **module_todos**
- **module_version**


---

<a href="../../nf_core/modules/lint/__init__.py#L32"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLintException`
Exception raised when there was an error with module linting 





---

<a href="../../nf_core/modules/lint/__init__.py#L38"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `LintResult`
An object to hold the results of a lint test 

<a href="../../nf_core/modules/lint/__init__.py#L41"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(mod, lint_test, message, file_path)
```









---

<a href="../../nf_core/modules/lint/__init__.py#L49"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleLint`
An object for linting modules either in a clone of the 'nf-core/modules' repository or in any nf-core pipeline directory 

<a href="../../nf_core/modules/lint/__init__.py#L65"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../nf_core/modules/lint/__init__.py#L264"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `filter_tests_by_key`

```python
filter_tests_by_key(key)
```

Filters the tests by the supplied key 

---

<a href="../../nf_core/modules/lint/__init__.py#L141"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_all_lint_tests`

```python
get_all_lint_tests(is_pipeline)
```





---

<a href="../../nf_core/modules/lint/__init__.py#L156"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint`

```python
lint(
    module=None,
    registry='quay.io',
    key=(),
    all_modules=False,
    print_results=True,
    show_passed=False,
    sort_by='test',
    local=False,
    fix_version=False
)
```

Lint all or one specific module 

First gets a list of all local modules (in modules/local/process) and all modules installed from nf-core (in modules/nf-core) 

For all nf-core modules, the correct file structure is assured and important file content is verified. If directory subject to linting is a clone of 'nf-core/modules', the files necessary for testing the modules are also inspected. 

For all local modules, the '.nf' file is checked for some important flags, and warnings are issued if untypical content is found. 

:param module:          A specific module to lint :param print_results:   Whether to print the linting results :param show_passed:     Whether passed tests should be shown as well :param fix_version:     Update the module version if a newer version is available :param hide_progress:   Don't show progress bars 

:returns:               A ModuleLint object containing information of  the passed, warned and failed tests 

---

<a href="../../nf_core/modules/lint/__init__.py#L308"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_module`

```python
lint_module(mod, progress_bar, registry, local=False, fix_version=False)
```

Perform linting on one module 

If the module is a local module we only check the `main.nf` file, and issue warnings instead of failures. 

If the module is a nf-core module we check for existence of the files 
- main.nf 
- meta.yml And verify that their content conform to the nf-core standards. 

If the linting is run for modules in the central nf-core/modules repo (repo_type==modules), files that are relevant for module testing are also examined 

---

<a href="../../nf_core/modules/lint/__init__.py#L279"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `lint_modules`

```python
lint_modules(modules, registry='quay.io', local=False, fix_version=False)
```

Lint a list of modules 



**Args:**
 
 - <b>`modules`</b> ([NFCoreModule]):  A list of module objects 
 - <b>`registry`</b> (str):  The container registry to use. Should be quay.io in most situations. 
 - <b>`local`</b> (boolean):  Whether the list consist of local or nf-core modules 
 - <b>`fix_version`</b> (boolean):  Fix the module version if a newer version is available 

---

<a href="../../nf_core/modules/lint/__init__.py#L455"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_summary`

```python
print_summary()
```

Print a summary table to the console. 

---

<a href="../../nf_core/modules/lint/__init__.py#L252"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `set_up_pipeline_files`

```python
set_up_pipeline_files()
```








---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
