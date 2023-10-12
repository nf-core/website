<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/create.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.create`

Creates a nf-core pipeline matching the current organization's specification based on a template.

---

<a href="../../../../../../tools/nf_core/create.py#L30"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `PipelineCreate`

Creates a nf-core pipeline a la carte from the nf-core best-practice template.

**Args:**

- <b>`name`</b> (str): Name for the pipeline.
- <b>`description`</b> (str): Description for the pipeline.
- <b>`author`</b> (str): Authors name of the pipeline.
- <b>`version`</b> (str): Version flag. Semantic versioning only. Defaults to `1.0dev`.
- <b>`no_git`</b> (bool): Prevents the creation of a local Git repository for the pipeline. Defaults to False.
- <b>`force`</b> (bool): Overwrites a given workflow directory with the same name. Defaults to False. May the force be with you.
- <b>`outdir`</b> (str): Path to the local output directory.

<a href="../../../../../../tools/nf_core/create.py#L44"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    name,
    description,
    author,
    version='1.0dev',
    no_git=False,
    force=False,
    outdir=None,
    template_yaml_path=None,
    plain=False
)
```

---

<a href="../../../../../../tools/nf_core/create.py#L89"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create_param_dict`

```python
create_param_dict(name, description, author, version, template_yaml_path, plain)
```

Creates a dictionary of parameters for the new pipeline.

**Args:**

- <b>`template_yaml_path`</b> (str): Path to YAML file containing template parameters.

---

<a href="../../../../../../tools/nf_core/create.py#L167"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `customize_template`

```python
customize_template(template_areas)
```

Customizes the template parameters.

**Args:**

- <b>`name`</b> (str): Name for the pipeline.
- <b>`description`</b> (str): Description for the pipeline.
- <b>`author`</b> (str): Authors name of the pipeline.

---

<a href="../../../../../../tools/nf_core/create.py#L478"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_pipeline_logo`

```python
download_pipeline_logo(url, img_fn)
```

Attempt to download a logo from the website. Retry if it fails.

---

<a href="../../../../../../tools/nf_core/create.py#L385"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `fix_linting`

```python
fix_linting()
```

Updates the .nf-core.yml with linting configurations for a customized pipeline.

---

<a href="../../../../../../tools/nf_core/create.py#L190"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_param`

```python
get_param(param_name, passed_value, template_yaml, template_yaml_path)
```

---

<a href="../../../../../../tools/nf_core/create.py#L519"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `git_init_pipeline`

```python
git_init_pipeline()
```

Initialises the new pipeline as a Git repository and submits first commit.

**Raises:**

- <b>`UserWarning`</b>: if Git default branch is set to 'dev' or 'TEMPLATE'.

---

<a href="../../../../../../tools/nf_core/create.py#L216"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `init_pipeline`

```python
init_pipeline()
```

Creates the nf-core pipeline.

---

<a href="../../../../../../tools/nf_core/create.py#L465"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_pipeline_logo`

```python
make_pipeline_logo()
```

Fetch a logo for the new pipeline from the nf-core website

---

<a href="../../../../../../tools/nf_core/create.py#L212"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_wf_author`

```python
prompt_wf_author()
```

---

<a href="../../../../../../tools/nf_core/create.py#L208"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_wf_description`

```python
prompt_wf_description()
```

---

<a href="../../../../../../tools/nf_core/create.py#L199"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `prompt_wf_name`

```python
prompt_wf_name()
```

---

<a href="../../../../../../tools/nf_core/create.py#L357"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `remove_nf_core_in_bug_report_template`

```python
remove_nf_core_in_bug_report_template()
```

Remove the field mentioning nf-core documentation in the github bug report template

---

<a href="../../../../../../tools/nf_core/create.py#L234"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `render_template`

```python
render_template()
```

Runs Jinja to create a new nf-core pipeline.

---

<a href="../../../../../../tools/nf_core/create.py#L334"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `update_nextflow_schema`

```python
update_nextflow_schema()
```

Removes unused parameters from the nextflow schema.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
