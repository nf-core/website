---
title: Deprecating components
subtitle: Deprecate a module or subworkflow
shortTitle: Deprecating components
---

## Deprecating a module or subworkflow

Sometimes modules or subworkflows become outdated and need to be deprecated (available, but no longer recommended).
They are not deleted straight away, as they could be used in private repositories, or used on other platforms.
Modules deprecated for more than one year are removed at the annual spring clean.

Once the alternative is available on nf-core modules, add a message explaining that the module/subworkflow is deprecated, and an `assert` in the code body to print it.

## Checklist

For modules:

- [ ] The replacement module exists in nf-core/modules
- [ ] `main.nf`: the deprecation message names the replacement as `nf-core/modules/<tool>/<subtool>`
- [ ] `main.nf`: `assert false: deprecation_message` comes first in `script:`, before any other check (e.g. the conda profile check)
- [ ] `main.nf`: the `stub:` block is removed
- [ ] `meta.yml`: `deprecated: true`, `deprecation_date` and `deprecation_pr` are set
- [ ] `tests/main.nf.test`: a single test with dummy inputs and no `setup` block asserts the failure and the deprecation message
- [ ] `tests/`: snapshots, test configs and any other test files are removed
- [ ] `nf-core modules lint <tool>/<subtool>` and the nf-test pass

For subworkflows:

- [ ] The replacement subworkflow exists in nf-core/subworkflows
- [ ] `main.nf`: the deprecation message names the replacement and `assert false: deprecation_message` comes first in `main:`
- [ ] `tests/main.nf.test`: the tests have no `setup` block and assert the failure

## Updating main.nf

For modules:

```nextflow title="main.nf"
process OLD_MODULE {
  ...

  script:
  def deprecation_message = """
  WARNING: This module has been deprecated. Please use nf-core/modules/<tool>/<subtool>

  Reason:
  This module is no longer fit for purpose because ...
  """
  assert false: deprecation_message
  ...
}
```

Remove the `stub:` block: with `-stub`, Nextflow runs `script:` for processes without a stub, so the `assert` also stops stub runs.

For subworkflows:

```nextflow title="main.nf"
...
  main:

  def deprecation_message = """
  WARNING: This subworkflow has been deprecated. Please use subworkflows/nf-core/<subworkflow>

  Reason:
  This subworkflow is no longer fit for purpose because ...
  """

  assert false: deprecation_message
}
```

The purpose of the `assert` is to introduce a mechanism which stops the pipeline and alerts the developer when
an automatic update of the module/subworkflow is performed.

## Updating meta.yml

Add `deprecated: true{:yml}` above `keywords:` in the module's `meta.yml`, with the date (`YYYY-MM-DD`) and the pull request of the deprecation:

```yaml title="meta.yml" {3-5}
name: old_module
description: ...
deprecated: true
deprecation_date: "2026-11-19"
deprecation_pr: "https://github.com/nf-core/modules/pull/1234"
keywords:
  - ...
```

`deprecated` is used by nf-core tools and the website to flag deprecated modules.
`deprecation_date` tells if the module can be removed at the spring clean.

## Updating the nf-tests

The `assert false` makes every test fail, so the tests only check that failure.
For modules, keep a single test with dummy inputs, without a `setup` block, and remove the snapshots and test configs:

```nextflow title="main.nf.test"
...
  test("deprecated") {

    when {
      process {
        """
        input[0] = [ [ id:'test' ], file('input') ]
        """
      }
    }

    then {
      assertAll(
        { assert process.failed },
        { assert process.errorReport.contains("WARNING: This module has been deprecated. Please use nf-core/modules/<tool>/<subtool>") }
      )
    }
  }
...
```

For subworkflows, remove anything in the `setup` block and update the `then` block:

```nextflow title="main.nf.test"
...
  then {
    assertAll(
      { assert workflow.failed }
    )
  }
...
```
