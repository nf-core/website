---
title: Deprecate modules
subtitle: How to deprecate a module or subworkflow
---

## Deprecate a module

Modules or subworkflows sometimes become outdated and need to be deprecated (available, but no longer recommended).
Don't delete deprecated modules or subworkflows. They may be used in private repositories or on other platforms.

Once an alternative is available in nf-core modules, follow this procedure:

1. Add a message to the top of the module code stating the module is deprecated
2. Add an `assert` statement in the code body to print a deprecation message

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

The `assert` statement stops the pipeline and alerts developers when they automatically update the module or subworkflow.

### Update the nf-tests

Adding `assert false` means the nf-test will always fail. Update the `then` block of the nf-test with the following code:

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
