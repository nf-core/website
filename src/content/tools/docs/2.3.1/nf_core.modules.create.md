<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/create.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.create`

The ModuleCreate class handles generating of module templates

---

<a href="../../../../../../tools/nf_core/modules/create.py#L27"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModuleCreate`

<a href="../../../../../../tools/nf_core/modules/create.py#L28"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    directory='.',
    tool='',
    author=None,
    process_label=None,
    has_meta=None,
    force=False,
    conda_name=None,
    conda_version=None,
    repo_type=None
)
```

---

<a href="../../../../../../tools/nf_core/modules/create.py#L60"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create`

```python
create()
```

Create a new DSL2 module from the nf-core template.

Tool should be named just <tool> or <tool/subtool> e.g fastqc or samtools/sort, respectively.

If <directory> is a pipeline, this function creates a file called: '<directory>/modules/local/tool.nf' OR '<directory>/modules/local/tool_subtool.nf'

If <directory> is a clone of nf-core/modules, it creates or modifies the following files:

modules/modules/tool/subtool/ _ main.nf _ meta.yml modules/tests/modules/tool/subtool/ _ main.nf _ test.yml \* nextflow.config tests/config/pytest_modules.yml

The function will attempt to automatically find a Bioconda package called <tool> and matching Docker / Singularity images from BioContainers.

---

<a href="../../../../../../tools/nf_core/modules/create.py#L294"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_module_dirs`

```python
get_module_dirs()
```

Given a directory and a tool/subtool, set the file paths and check if they already exist

Returns dict: keys are relative paths to template files, vals are target paths.

---

<a href="../../../../../../tools/nf_core/modules/create.py#L271"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `render_template`

```python
render_template()
```

Create new module files with Jinja2.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
