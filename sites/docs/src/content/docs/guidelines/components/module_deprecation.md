---
title: Deprecating a module
subtitle: How to deprecate a module or subworkflow
---

## Deprecating a module

Sometimes modules or subworkflows become outdated and need to be deprecated (available, but no longer recommended).
These modules or subworkflows should not be deleted as they could be used in private repositories, or used on other
platforms. The recommended procedure is, once the alternative is available on nf-core modules, add a message to the
top of the module code saying this module is deprecated, and an `assert` in the code body to print a deprecation
message like so:

```groovy title="main.nf"
def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

Reason:
This module is no longer fit for purpose because ...

"""

process OLD_MODULE {
  ...

  script:
  assert false: deprecation_message
}
```

The purpose of the `assert` is to introduce a mechanism which stops the pipeline and alerts the developer when
an automatic update of the module/subworkflow is performed.
