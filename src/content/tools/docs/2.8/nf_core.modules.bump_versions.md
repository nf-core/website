<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.bump_versions`

Bump versions for all modules on nf-core/modules or for a single module

---

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L28"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleVersionBumper`

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L29"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline_dir, remote_url=None, branch=None, no_pull=False)
```

---

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L117"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `bump_module_version`

```python
bump_module_version(module)
```

Bump the bioconda and container version of a single NFCoreModule

**Args:**

- <b>`module`</b>: NFCoreModule

---

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L38"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `bump_versions`

```python
bump_versions(module=None, all_modules=False, show_uptodate=False)
```

Bump the container and conda version of single module or all modules

Looks for a bioconda tool version in the `main.nf` file of the module and checks whether are more recent version is available. If yes, then tries to get docker/singularity container links and replace the bioconda version and the container links in the main.nf file of the respective module.

**Args:**

- <b>`module`</b>: a specific module to update
- <b>`all_modules`</b>: whether to bump versions for all modules

---

<a href="../../../../../../tools/nf_core/modules/bump_versions.py#L226"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_bioconda_version`

```python
get_bioconda_version(module)
```

Extract the bioconda version from a module

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
