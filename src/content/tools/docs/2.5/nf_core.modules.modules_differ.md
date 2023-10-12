<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_differ`

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesDiffer`

Static class that provides functionality for computing diffs between different instances of a module

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L186"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `append_modules_json_diff`

```python
append_modules_json_diff(
    diff_path,
    old_modules_json,
    new_modules_json,
    modules_json_path,
    for_git=True
)
```

Compare the new modules.json and builds a diff

**Args:**

- <b>`diff_fn`</b> (str): The diff file to be appended
- <b>`old_modules_json`</b> (nested dict): The old modules.json
- <b>`new_modules_json`</b> (nested dict): The new modules.json
- <b>`modules_json_path`</b> (str): The path to the modules.json
- <b>`for_git`</b> (bool): indicates whether the diff file is to be compatible with `git apply`. If true it adds a/ and b/ prefixes to the file paths

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_diffs`

```python
get_module_diffs(
    from_dir,
    to_dir,
    for_git=True,
    dsp_from_dir=None,
    dsp_to_dir=None
)
```

Compute the diff between the current module version and the new version.

**Args:**

- <b>`from_dir`</b> (strOrPath): The folder containing the old module files
- <b>`to_dir`</b> (strOrPath): The folder containing the new module files
- <b>`path_in_diff`</b> (strOrPath): The directory displayed containing the module file in the diff. Added so that temporary dirs are not shown
- <b>`for_git`</b> (bool): indicates whether the diff file is to be compatible with `git apply`. If true it adds a/ and b/ prefixes to the file paths
- <b>`dsp_from_dir`</b> (str | Path): The from directory to display in the diff
- <b>`dsp_to_dir`</b> (str | Path): The to directory to display in the diff

**Returns:**

- <b>`dict[str, (ModulesDiffer.DiffEnum, str)]`</b>: A dictionary containing the diff type and the diff string (empty if no diff)

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L321"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_new_and_old_lines`

```python
get_new_and_old_lines(patch)
```

Parse a patch for a file, and return the contents of the modified parts for both the old and new versions

**Args:**

- <b>`patch`</b> (str): The patch in unified diff format

**Returns:**

- <b>`([[str]], [[str]])`</b>: Lists of old and new lines for each hunk (modified part the file)

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L270"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `per_file_patch`

```python
per_file_patch(patch_fn)
```

Splits a patch file for several files into one patch per file.

**Args:**

- <b>`patch_fn`</b> (str | Path): The path to the patch file

**Returns:**

- <b>`dict[str, str]`</b>: A dictionary indexed by the filenames with the file patches as values

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L220"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_diff`

```python
print_diff(
    module,
    repo_name,
    from_dir,
    to_dir,
    current_version=None,
    new_version=None,
    dsp_from_dir=None,
    dsp_to_dir=None
)
```

Prints the diffs between two module versions to the terminal

**Args:**

- <b>`module`</b> (str): The module name
- <b>`repo_name`</b> (str): The name of the repo where the module resides
- <b>`from_dir`</b> (str | Path): The directory containing the old module files
- <b>`to_dir`</b> (str | Path): The directory containing the new module files
- <b>`module_dir`</b> (str): The path to the current installation of the module
- <b>`current_version`</b> (str): The installed version of the module
- <b>`new_version`</b> (str): The version of the module the diff is computed against
- <b>`dsp_from_dir`</b> (str | Path): The 'from' directory displayed in the diff
- <b>`dsp_to_dir`</b> (str | Path): The 'to' directory displayed in the diff

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L425"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `try_apply_patch`

```python
try_apply_patch(module, repo_name, patch_path, module_dir, reverse=False)
```

Try applying a full patch file to a module

**Args:**

- <b>`module`</b> (str): Name of the module
- <b>`repo_name`</b> (str): Name of the repository where the module resides
- <b>`patch_path`</b> (str): The absolute path to the patch file to be applied
- <b>`module_dir`</b> (Path): The directory containing the module

**Returns:**

- <b>`dict[str, str]`</b>: A dictionary with file paths (relative to the pipeline dir) as keys and the patched file contents as values

**Raises:**

- <b>`LookupError`</b>: If the patch application fails in a file

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L362"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `try_apply_single_patch`

```python
try_apply_single_patch(file_lines, patch, reverse=False)
```

Tries to apply a patch to a modified file. Since the line numbers in the patch does not agree if the file is modified, the old and new lines in the patch are reconstructed and then we look for the old lines in the modified file. If all hunk in the patch are found in the new file it is updated with the new lines from the patch file.

**Args:**

- <b>`new_fn`</b> (str | Path): Path to the modified file
- <b>`patch`</b> (str | Path): (Outdated) patch for the file
- <b>`reverse`</b> (bool): Apply the patch in reverse

**Returns:**

- <b>`[str]`</b>: The patched lines of the file

**Raises:**

- <b>`LookupError`</b>: If it fails to find the old lines from the patch in the file.

---

<a href="../../../../../../tools/nf_core/modules/modules_differ.py#L123"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `write_diff_file`

```python
write_diff_file(
    diff_path,
    module,
    repo_name,
    from_dir,
    to_dir,
    current_version=None,
    new_version=None,
    file_action='a',
    for_git=True,
    dsp_from_dir=None,
    dsp_to_dir=None
)
```

Writes the diffs of a module to the diff file.

**Args:**

- <b>`diff_path`</b> (str | Path): The path to the file that should be appended
- <b>`module`</b> (str): The module name
- <b>`repo_name`</b> (str): The name of the repo where the module resides
- <b>`from_dir`</b> (str | Path): The directory containing the old module files
- <b>`to_dir`</b> (str | Path): The directory containing the new module files
- <b>`diffs`</b> (dict[str, (ModulesDiffer.DiffEnum, str)]): A dictionary containing the type of change and the diff (if any)
- <b>`module_dir`</b> (str | Path): The path to the current installation of the module
- <b>`current_version`</b> (str): The installed version of the module
- <b>`new_version`</b> (str): The version of the module the diff is computed against
- <b>`for_git`</b> (bool): indicates whether the diff file is to be compatible with `git apply`. If true it adds a/ and b/ prefixes to the file paths
- <b>`dsp_from_dir`</b> (str | Path): The 'from' directory displayed in the diff
- <b>`dsp_to_dir`</b> (str | Path): The 'to' directory displayed in the diff

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
