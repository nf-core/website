<!-- markdownlint-disable -->

<a href="../../nf_core/lint_utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint_utils`





---

<a href="../../nf_core/lint_utils.py#L19"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `print_joint_summary`

```python
print_joint_summary(lint_obj, module_lint_obj)
```

Print a joint summary of the general pipe lint tests and the module lint tests 


---

<a href="../../nf_core/lint_utils.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `print_fixes`

```python
print_fixes(lint_obj)
```

Prints available and applied fixes 


---

<a href="../../nf_core/lint_utils.py#L57"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `run_prettier_on_file`

```python
run_prettier_on_file(file)
```

Run the pre-commit hook prettier on a file. 



**Args:**
 
 - <b>`file`</b> (Path | str):  A file identifier as a string or pathlib.Path. 

Warns: If Prettier is not installed, a warning is logged. 


---

<a href="../../nf_core/lint_utils.py#L88"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `dump_json_with_prettier`

```python
dump_json_with_prettier(file_name, file_content)
```

Dump a JSON file and run prettier on it. 

**Args:**
 
 - <b>`file_name`</b> (Path | str):  A file identifier as a string or pathlib.Path. 
 - <b>`file_content`</b> (dict):  Content to dump into the JSON file 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
