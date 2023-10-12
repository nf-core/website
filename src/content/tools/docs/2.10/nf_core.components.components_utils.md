<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/components/components_utils.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.components.components_utils`

---

<a href="../../../../../../tools/nf_core/components/components_utils.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_repo_info`

```python
get_repo_info(directory, use_prompt=True)
```

Determine whether this is a pipeline repository or a clone of nf-core/modules

---

<a href="../../../../../../tools/nf_core/components/components_utils.py#L83"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `prompt_component_version_sha`

```python
prompt_component_version_sha(
    component_name,
    component_type,
    modules_repo,
    installed_sha=None
)
```

Creates an interactive questionary prompt for selecting the module/subworkflow version

**Args:**

- <b>`component_name`</b> (str): Module/subworkflow name,
- <b>`component_type`</b> (str): "modules" or "subworkflows",
- <b>`modules_repo`</b> (ModulesRepo): Modules repo the module/subworkflow originate in
- <b>`installed_sha`</b> (str): Optional extra argument to highlight the current installed version

**Returns:**

- <b>`git_sha`</b> (str): The selected version of the module/subworkflow

---

<a href="../../../../../../tools/nf_core/components/components_utils.py#L129"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_components_to_install`

```python
get_components_to_install(subworkflow_dir)
```

Parse the subworkflow main.nf file to retrieve all imported modules and subworkflows.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
