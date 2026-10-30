---
title: Deprecating components
subtitle: Deprecate a module or subworkflow
shortTitle: Deprecating components
---

## Deprecating a module or subworkflow

Sometimes modules or subworkflows become outdated and need to be deprecated (available, but no longer recommended).
These modules or subworkflows should not be deleted as they could be used in private repositories, or used on other
platforms.

The recommended procedure is, once the alternative is available on nf-core modules, to add a message explaining that
the module/subworkflow is deprecated, and an `assert` in the code body to print a deprecation
message like.

For modules:

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

For subworkflows:

```groovy title="main.nf"
...
  main:

  def deprecation_message = """
  WARNING: This subworkflow has been deprecated. Please use subworkflows/nf-core//path/to/new/module

  Reason:
  This subworkflow is no longer fit for purpose because ...
  """

  assert false: deprecation_message
}
```

The purpose of the `assert` is to introduce a mechanism which stops the pipeline and alerts the developer when
an automatic update of the module/subworkflow is performed.

### Updating meta.yml

Add `deprecated: true{:yml}` above `keywords:` to the module's `meta.yml` to mark it as deprecated in the nf-core registry and tooling:

```yaml title="meta.yml" {3}
name: old_module
description: ...
deprecated: true
keywords:
  - ...
```

This field is used by nf-core tools and the website to flag deprecated modules.

### Updating the nf-tests

Remove anything in the setup block.
The addition of `assert false` will mean the nf-test will always fail.
The `then` part of the nf-test should be updated with the following code for modules:

```groovy title="main.nf.test"
...
  then {
    assertAll(
      { assert process.failed },
      { assert process.errorReport.contains("WARNING: This module has been deprecated.")}
    )
  }
...
```

Or the following code for subworkflows:

```groovy title="main.nf.test"
...
  then {
    assertAll(
      { assert workflow.failed }
    )
  }
...
```
