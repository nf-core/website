<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/create.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.create`

Creates a nf-core pipeline matching the current organization's specification based on a template.

---

<a href="../../../../../../tools/nf_core/create.py#L23"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../../../../../../tools/nf_core/create.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    name,
    description,
    author,
    version='1.0dev',
    no_git=False,
    force=False,
    outdir=None
)
```

---

<a href="../../../../../../tools/nf_core/create.py#L161"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_pipeline_logo`

```python
download_pipeline_logo(url, img_fn)
```

Attempt to download a logo from the website. Retry if it fails.

---

<a href="../../../../../../tools/nf_core/create.py#L202"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `git_init_pipeline`

```python
git_init_pipeline()
```

Initialises the new pipeline as a Git repository and submits first commit.

---

<a href="../../../../../../tools/nf_core/create.py#L51"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `init_pipeline`

```python
init_pipeline()
```

Creates the nf-core pipeline.

---

<a href="../../../../../../tools/nf_core/create.py#L148"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `make_pipeline_logo`

```python
make_pipeline_logo()
```

Fetch a logo for the new pipeline from the nf-core website

---

<a href="../../../../../../tools/nf_core/create.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `render_template`

```python
render_template()
```

Runs Jinja to create a new nf-core pipeline.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
