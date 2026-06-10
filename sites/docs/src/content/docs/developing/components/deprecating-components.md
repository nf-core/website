---
title: Deprecating modules
subtitle: Deprecate a module or subworkflow
shortTitle: Deprecating modules
---

## Deprecating a module

Sometimes modules or subworkflows become outdated and need to be deprecated (available, but no longer recommended).
These modules or subworkflows should not be deleted as they could be used in private repositories, or used on other
platforms.
The recommended procedure is, once the alternative is available on nf-core modules, add a message to the
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
The `then` part of the nf-test should be updated with the following code:

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
