<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/module_test.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.module_test`

The ModulesTest class runs the tests locally

---

<a href="../../../../../../tools/nf_core/modules/module_test.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesTest`

Class to run module pytests.

...

Attributes
---------- module_name : str name of the tool to run tests for no_prompts : bool flat indicating if prompts are used pytest_args : tuple additional arguments passed to pytest command

Methods
------- run(): Run test steps \_check_inputs(): Check inputs. Ask for module_name if not provided and check that the directory exists \_set_profile(): Set software profile \_run_pytests(self): Run pytest

<a href="../../../../../../tools/nf_core/modules/module_test.py#L51"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    module_name=None,
    no_prompts=False,
    pytest_args='',
    remote_url=None,
    branch=None,
    no_pull=False
)
```

---

<a href="../../../../../../tools/nf_core/modules/module_test.py#L66"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `run`

```python
run()
```

Run test steps

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
