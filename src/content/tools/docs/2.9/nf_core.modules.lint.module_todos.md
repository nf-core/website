<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/modules/lint/module_todos.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.lint.module_todos`

---

<a href="../../../../../../tools/nf_core/modules/lint/module_todos.py#L8"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `module_todos`

```python
module_todos(_, module)
```

Look for TODO statements in the module files

The nf-core module template contains a number of comment lines to help developers of new modules know where they need to edit files and add content. They typically have the following format:

.. code-block:: groovy

// TODO nf-core: Make some kind of change to the workflow here

..or in markdown:

.. code-block:: html

 <!-- TODO nf-core: Add some detail to the docs here -->

This lint test runs through all files in the module and searches for these lines. If any are found they will throw a warning.

.. tip:: Note that many GUI code editors have plugins to list all instances of _TODO_ in a given project directory. This is a very quick and convenient way to get started on your pipeline!

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
