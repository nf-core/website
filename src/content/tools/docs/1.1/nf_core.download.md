<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.download`

Download a nf-core pipeline

---

<a href="../../../../../../tools/nf_core/download.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DownloadWorkflow`

<a href="../../../../../../tools/nf_core/download.py#L23"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(pipeline, release=None, singularity=False, outdir=None)
```

Set class variables

---

<a href="../../../../../../tools/nf_core/download.py#L171"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_shub_image`

```python
download_shub_image(container)
```

Download singularity images from singularity-hub

---

<a href="../../../../../../tools/nf_core/download.py#L147"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_wf_files`

```python
download_wf_files()
```

Download workflow files from GitHub - save in outdir

---

<a href="../../../../../../tools/nf_core/download.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `download_workflow`

```python
download_workflow()
```

Main function to download a nf-core workflow

---

<a href="../../../../../../tools/nf_core/download.py#L79"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `fetch_workflow_details`

```python
fetch_workflow_details(wfs)
```

Fetch details of nf-core workflow to download

params:

- wfs A nf_core.list.Workflows object

---

<a href="../../../../../../tools/nf_core/download.py#L160"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `find_singularity_images`

```python
find_singularity_images()
```

Find singularity image names for workflow

---

<a href="../../../../../../tools/nf_core/download.py#L217"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `pull_singularity_image`

```python
pull_singularity_image(container)
```

Use a local installation of singularity to pull an image from docker hub

---

<a href="../../../../../../tools/nf_core/download.py#L237"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate_md5`

```python
validate_md5(fname, expected)
```

Calculate the md5sum for a file on the disk and validate with expected

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
