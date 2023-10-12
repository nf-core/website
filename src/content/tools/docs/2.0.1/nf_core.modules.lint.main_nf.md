<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/main_nf.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.main_nf`

Lint the main.nf file of a module

---

<a href="../../../../../../tools/nf_core/modules/lint/main_nf.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main_nf`

```python
main_nf(module_lint_object, module)
```

Lint a single main.nf module file Can also be used to lint local module files, in which case failures should be interpreted as warnings

---

<a href="../../../../../../tools/nf_core/modules/lint/main_nf.py#L116"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_script_section`

```python
check_script_section(lines)
```

Lint the script section Checks whether 'def sotware' and 'def prefix' are defined

---

<a href="../../../../../../tools/nf_core/modules/lint/main_nf.py#L137"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_process_section`

```python
check_process_section(lines)
```

Lint the section of a module between the process definition and the 'input:' definition Specifically checks for correct software versions and containers

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
