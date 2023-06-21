<!-- markdownlint-disable -->

<a href="../../nf_core/lint/pipeline_todos.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.pipeline_todos`





---

<a href="../../nf_core/lint/pipeline_todos.py#L9"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `pipeline_todos`

```python
pipeline_todos(root_dir=None)
```

Check for nf-core *TODO* lines. 

The nf-core workflow template contains a number of comment lines to help developers of new pipelines know where they need to edit files and add content. They typically have the following format: 

.. code-block:: groovy 

 // TODO nf-core: Make some kind of change to the workflow here 

..or in markdown: 

.. code-block:: html 

 <!-- TODO nf-core: Add some detail to the docs here --> 

This lint test runs through all files in the pipeline and searches for these lines. If any are found they will throw a warning. 

.. tip:: Note that many GUI code editors have plugins to list all instances of *TODO*  in a given project directory. This is a very quick and convenient way to get  started on your pipeline! 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
