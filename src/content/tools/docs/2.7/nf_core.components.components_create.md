<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/components/components_create.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.components.components_create`

---

<a href="../../../../../../tools/nf_core/components/components_create.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `render_template`

```python
render_template(component_type, object_attrs, file_paths)
```

Create new module/subworkflow files with Jinja2.

---

<a href="../../../../../../tools/nf_core/components/components_create.py#L43"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `collect_name_prompt`

```python
collect_name_prompt(name, component_type)
```

Collect module/subworkflow info via prompt if empty or invalid

---

<a href="../../../../../../tools/nf_core/components/components_create.py#L91"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_component_dirs`

```python
get_component_dirs(
    component_type,
    repo_type,
    directory,
    org,
    name,
    supername,
    subname,
    new_dir,
    force_overwrite
)
```

Given a directory and a tool/subtool or subworkflow, set the file paths and check if they already exist

Returns dict: keys are relative paths to template files, vals are target paths.

---

<a href="../../../../../../tools/nf_core/components/components_create.py#L154"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_username`

```python
get_username(author)
```

Prompt for GitHub username

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
