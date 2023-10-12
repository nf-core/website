<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/mulled.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.mulled`

Generate the name of a BioContainers mulled image version 2.

---

<a href="../../../../../../tools/nf_core/modules/mulled.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `MulledImageNameGenerator`

Define a service class for generating BioContainers version 2 mulled image names.

Adapted from https://gist.github.com/natefoo/19cefeedd1942c30f9d88027a61b3f83.

---

<a href="../../../../../../tools/nf_core/modules/mulled.py#L50"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `generate_image_name`

```python
generate_image_name(
    targets: Iterable[Tuple[str, str]],
    build_number: int = 0
) → str
```

Generate the name of a BioContainers mulled image version 2.

**Args:**

- <b>`targets`</b>: One or more tool, version pairs of the multi-tool container image.
- <b>`build_number`</b>: The build number for this image. This is an incremental value that starts at zero.

---

<a href="../../../../../../tools/nf_core/modules/mulled.py#L62"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `image_exists`

```python
image_exists(image_name: str) → bool
```

Check whether a given BioContainers image name exists via a call to the quay.io API.

---

<a href="../../../../../../tools/nf_core/modules/mulled.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `parse_targets`

```python
parse_targets(specifications: Iterable[str]) → List[Tuple[str, str]]
```

Parse tool, version pairs from specification strings.

**Args:**

- <b>`specifications`</b>: An iterable of strings that contain tools and their versions.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
