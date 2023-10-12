<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.workflow.parameters`

## **Global Variables**

- **NFCORE_PARAMS_SCHEMA_URI**

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Parameters`

Contains a static factory method for :class:`Parameter` object creation.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create_from_json`

```python
create_from_json(parameters_json, schema_json='')
```

Creates a list of Parameter objects from a description in JSON.

**Args:**

- <b>`parameters_json`</b> (str): Parameter(s) description in JSON.
- <b>`schema`</b> (str): Parameter schema in JSON.

**Returns:**

- <b>`list`</b>: Parameter objects.

**Raises:**

- <b>`ValidationError`</b>: When the parameter JSON violates the schema.
- <b>`LookupError`</b>: When the schema cannot be downloaded.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L80"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `in_full_json`

```python
in_full_json(parameters, indent=0)
```

Converts a list of Parameter objects into JSON. All attributes are written.

**Args:**

- <b>`parameters`</b> (list): List of :class:`Parameter` objects.
- <b>`indent`</b> (integer): String output indentation. Defaults to 0.

**Returns:**

- <b>`list`</b>: JSON formatted parameters.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `in_nextflow_json`

```python
in_nextflow_json(parameters, indent=0)
```

Converts a list of Parameter objects into JSON, readable by Nextflow.

**Args:**

- <b>`parameters`</b> (list): List of :class:`Parameter` objects.
- <b>`indent`</b> (integer): String output indentation. Defaults to 0.

**Returns:**

- <b>`list`</b>: JSON formatted parameters.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L106"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Parameter`

Holds information about a workflow parameter.

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L110"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(param_builder)
```

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L131"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `as_dict`

```python
as_dict()
```

Describes its attibutes in a dictionary.

**Returns:**

- <b>`dict`</b>: Parameter object as key value pairs.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `builder`

```python
builder()
```

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L145"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

Validates the parameter's value. If the value is within the parameter requirements, no exception is thrown.

**Raises:**

- <b>`LookupError`</b>: Raised when no matching validator can be determined.
- <b>`AttributeError`</b>: Raised with description, if a parameter value violates the parameter constrains.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ParameterBuilder`

Parameter builder.

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L162"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L257"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `arity`

```python
arity(arity)
```

Sets the parameter regex pattern.

**Args:**

- <b>`pattern`</b> (str): Parameter regex pattern.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L280"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build`

```python
build()
```

Builds parameter object.

**Returns:**

- <b>`Parameter`</b>: Fresh from the factory.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L221"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `choices`

```python
choices(choices)
```

Sets the parameter value choices.

**Args:**

- <b>`choices`</b> (list): Parameter value choices.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L239"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `default`

```python
default(default)
```

Sets the parameter default value.

**Args:**

- <b>`default`</b> (str): Parameter default value.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L176"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `group`

```python
group(group)
```

Sets the parameter group tag

**Args:**

- <b>`group`</b> (str): Parameter group tag.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L194"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `label`

```python
label(label)
```

Sets the parameter label.

**Args:**

- <b>`label`</b> (str): Parameter label.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L185"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `name`

```python
name(name)
```

Sets the parameter name.

**Args:**

- <b>`name`</b> (str): Parameter name.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L230"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `param_type`

```python
param_type(param_type)
```

Sets the parameter type.

**Args:**

- <b>`param_type`</b> (str): Parameter type.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L248"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `pattern`

```python
pattern(pattern)
```

Sets the parameter regex pattern.

**Args:**

- <b>`pattern`</b> (str): Parameter regex pattern.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L266"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `render`

```python
render(render)
```

Sets the parameter render type.

**Args:**

- <b>`render`</b> (str): UI render type.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L275"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `required`

```python
required(required)
```

Sets the required parameter flag.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L203"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `usage`

```python
usage(usage)
```

Sets the parameter usage.

**Args:**

- <b>`usage`</b> (str): Parameter usage description.

---

<a href="../../../../../../tools/nf_core/workflow/parameters.py#L212"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `value`

```python
value(value)
```

Sets the parameter value.

**Args:**

- <b>`value`</b> (str): Parameter value.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
