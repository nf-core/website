<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_repo`

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesRepo`

An object to store details about the repository being used for modules.

Used by the `nf-core modules` top-level command with -r and -b flags, so that this can be used in the same way by all sub-commands.

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L19"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(repo='nf-core/modules', branch='master')
```

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_gh_file`

```python
download_gh_file(dl_filename, api_url)
```

Download a file from GitHub using the GitHub API

**Args:**

- <b>`dl_filename`</b> (string): Path to save file to
- <b>`api_url`</b> (string): GitHub API URL for file

**Raises:**
If a problem, raises an error

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L92"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_file_urls`

```python
get_module_file_urls(module, commit='')
```

Fetch list of URLs for a specific module

Takes the name of a module and iterates over the GitHub repo file tree. Loops over items that are prefixed with the path 'modules/<module_name>' and ignores anything that's not a blob. Also ignores the test/ subfolder.

Returns a dictionary with keys as filenames and values as GitHub API URLs. These can be used to then download file contents.

**Args:**

- <b>`module`</b> (string): Name of module for which to fetch a set of URLs

**Returns:**

- <b>`dict`</b>: Set of files and associated URLs as follows:

{

- <b>`'modules/fastqc/main.nf'`</b>: 'https://api.github.com/repos/nf-core/modules/git/blobs/65ba598119206a2b851b86a9b5880b5476e263c3',
- <b>`'modules/fastqc/meta.yml'`</b>: 'https://api.github.com/repos/nf-core/modules/git/blobs/0d5afc23ba44d44a805c35902febc0a382b17651' }

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_modules_file_tree`

```python
get_modules_file_tree()
```

Fetch the file list from the repo, using the GitHub API

Sets self.modules_file_tree self.modules_current_hash self.modules_avail_module_names

---

<a href="../../../../../../tools/nf_core/modules/modules_repo.py#L35"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `verify_modules_repo`

```python
verify_modules_repo()
```

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
