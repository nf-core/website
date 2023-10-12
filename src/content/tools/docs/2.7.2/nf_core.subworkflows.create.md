<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/subworkflows/create.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.create`

The SubworkflowCreate class handles generating of subworkflow templates

---

<a href="../../../../../../tools/nf_core/subworkflows/create.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SubworkflowCreate`

<a href="../../../../../../tools/nf_core/subworkflows/create.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(directory='.', subworkflow='', author=None, force=False)
```

---

<a href="../../../../../../tools/nf_core/subworkflows/create.py#L35"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create`

```python
create()
```

Create a new subworkflow from the nf-core template.

The subworkflow should be named as the main file type it operates on and a short description of the task performed e.g bam_sort or bam_sort_samtools, respectively.

If <directory> is a pipeline, this function creates a file called: '<directory>/subworkflows/local/subworkflow_name.nf'

If <directory> is a clone of nf-core/modules, it creates or modifies the following files:

subworkflows/nf-core/subworkflow*name/ * main.nf _ meta.yml tests/subworkflows/nf-core/subworkflow_name/ _ main.nf \_ test.yml \* nextflow.config tests/config/pytest_modules.yml

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
