<!-- markdownlint-disable -->

<a href="../../nf_core/modules/lint/main_nf.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.main_nf`
Lint the main.nf file of a module 


---

<a href="../../nf_core/modules/lint/main_nf.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main_nf`

```python
main_nf(module_lint_object, module, fix_version, registry, progress_bar)
```

Lint a ``main.nf`` module file 

Can also be used to lint local module files, in which case failures will be reported as warnings. 

The test checks for the following: 

* Software versions and containers are valid * The module has a process label and it is among  the standard ones. * If a ``meta`` map is defined as one of the modules  inputs it should be defined as one of the outputs,  and be correctly configured in the ``saveAs`` function. * The module script section should contain definitions  of ``software`` and ``prefix`` 


---

<a href="../../nf_core/modules/lint/main_nf.py#L171"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_script_section`

```python
check_script_section(lines)
```

Lint the script section Checks whether `def prefix` is defined and whether getProcessName is used for `versions.yml`. 


---

<a href="../../nf_core/modules/lint/main_nf.py#L192"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_when_section`

```python
check_when_section(lines)
```

Lint the when: section Checks whether the line is modified from 'task.ext.when == null || task.ext.when' 


---

<a href="../../nf_core/modules/lint/main_nf.py#L212"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_process_section`

```python
check_process_section(lines, registry, fix_version, progress_bar)
```

Lint the section of a module between the process definition and the 'input:' definition Specifically checks for correct software versions and containers 



**Args:**
 
 - <b>`lines`</b> (List[str]):  Content of process. 
 - <b>`registry`</b> (str):  Base Docker registry for containers. Typically quay.io. 
 - <b>`fix_version`</b> (bool):  Fix software version 
 - <b>`progress_bar`</b> (ProgressBar):  Progress bar to update. 



**Returns:**
 
 - <b>`Optional[bool]`</b>:  True if singularity and docker containers match, False otherwise. If process definition does not exist, None. 


---

<a href="../../nf_core/modules/lint/main_nf.py#L436"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_process_labels`

```python
check_process_labels(lines)
```








---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
