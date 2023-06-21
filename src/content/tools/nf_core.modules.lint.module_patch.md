<!-- markdownlint-disable -->

<a href="../../nf_core/modules/lint/module_patch.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.module_patch`





---

<a href="../../nf_core/modules/lint/module_patch.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `module_patch`

```python
module_patch(module_lint_obj, module: NFCoreModule)
```

Lint a patch file found in a module 

Checks that the file name is well formed, and that the patch can be applied in reverse with the correct result. 


---

<a href="../../nf_core/modules/lint/module_patch.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_patch_valid`

```python
check_patch_valid(module, patch_path)
```

Checks whether a patch is valid. Looks for lines like 
--- <path> +++ <path> @@ n,n n,n @@ and make sure that the come in the right order and that the reported paths exists. If the patch file performs file creation or deletion we issue a lint warning. 



**Args:**
 
 - <b>`module`</b> (NFCoreModule):  The module currently being linted 
 - <b>`patch_path`</b> (Path):  The absolute path to the patch file. 



**Returns:**
 
 - <b>`(bool)`</b>:  False if any test failed, True otherwise 


---

<a href="../../nf_core/modules/lint/module_patch.py#L152"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `patch_reversible`

```python
patch_reversible(module_lint_object, module, patch_path)
```

Try applying a patch in reverse to see if it is up to date 



**Args:**
 
 - <b>`module`</b> (NFCoreModule):  The module currently being linted 
 - <b>`patch_path`</b> (Path):  The absolute path to the patch file. 



**Returns:**
 
 - <b>`(bool)`</b>:  False if any test failed, True otherwise 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
