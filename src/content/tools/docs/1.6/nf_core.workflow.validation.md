<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/workflow/validation.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.workflow.validation`

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L12"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Validators`

Gives access to a factory method for objects of instance :class:`Validator` which returns the correct Validator for a given parameter type.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L20"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_validator_for_param`

```python
get_validator_for_param(parameter)
```

Determines matching :class:`Validator` instance for a given parameter.

**Returns:**

- <b>`Validator`</b>: Matching validator for a given :class:`Parameter`.

**Raises:**

- <b>`LookupError`</b>: In case no matching validator for a given parameter type can be determined.

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L43"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Validator`

Abstract base class for different parameter validators.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(parameter)
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L52"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L57"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `IntegerValidator`

Implementation for parameters of type integer.

**Args:**

- <b>`parameter`</b> (:class:`Parameter`): A Parameter object.

**Raises:**

- <b>`AttributeError`</b>: In case the argument is not of instance integer.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L67"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(parameter)
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L70"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

Validates an parameter integer value against a given range (choices). If the value is valid, no error is risen.

**Raises:**

- <b>`AtrributeError`</b>: Description of the value error.

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L90"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `StringValidator`

Implementation for parameters of type string.

**Args:**

- <b>`parameter`</b> (:class:`Parameter`): A Parameter object.

**Raises:**

- <b>`AttributeError`</b>: In case the argument is not of instance string.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L100"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(parameter)
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L103"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

Validates an parameter integer value against a given range (choices). If the value is valid, no error is risen.

**Raises:**

- <b>`AtrributeError`</b>: Description of the value error.

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L133"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `BooleanValidator`

Implementation for parameters of type boolean.

**Args:**

- <b>`parameter`</b> (:class:`Parameter`): A Parameter object.

**Raises:**

- <b>`AttributeError`</b>: In case the argument is not of instance boolean.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L143"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(parameter)
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L146"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

Validates an parameter boolean value. If the value is valid, no error is risen.

**Raises:**

- <b>`AtrributeError`</b>: Description of the value error.

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L159"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DecimalValidator`

Implementation for parameters of type boolean.

**Args:**

- <b>`parameter`</b> (:class:`Parameter`): A Parameter object.

**Raises:**

- <b>`AttributeError`</b>: In case the argument is not of instance decimal.

<a href="../../../../../../tools/nf_core/workflow/validation.py#L169"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(parameter)
```

---

<a href="../../../../../../tools/nf_core/workflow/validation.py#L172"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `validate`

```python
validate()
```

Validates an parameter boolean value. If the value is valid, no error is risen.

**Raises:**

- <b>`AtrributeError`</b>: Description of the value error.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
