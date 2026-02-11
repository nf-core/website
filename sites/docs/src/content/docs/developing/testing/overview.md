---
title: Overview
subtitle: Testing nf-core components with nf-test
shortTitle: Overview
weight: 1
---

nf-test provides a testing framework for Nextflow pipelines, modules, and workflows.
Snapshots form the foundation of nf-test, comparing current outputs against reference files to ensure your components work correctly and prevent regressions.

The [nf-test documentation](https://code.askimed.com/nf-test/docs/getting-started/) provides comprehensive information about the framework.

## Snapshots

Snapshots are used to compare the current output of a process, workflow, or function against a reference snapshot file (`*.nf.test.snap`).

### Using snapshots

Create snapshots using the `snapshot` keyword. The `match` method checks if the snapshot corresponds to the expected data in the snap file:

```groovy
// Create a snapshot of a workflow channel
assert snapshot(workflow.out.channel1).match('channel1')

// Snapshot all output channels of a process
assert snapshot(process.out).match()

// Snapshot a specific file
assert snapshot(path(process.out.get(0))).match()

// Snapshot the result of a function
assert snapshot(function.result).match()
```

The first test run generates a JSON snapshot file. Subsequent runs compare against this file. Commit snapshot files with code changes and review them during code review.

### Assigning parameters or configs

You can specify parameters or include config files in your tests.

```groovy
params {
  foo = 'bar'
}
```

```groovy
config "./nextflow.config"
```

Use `withName` selectors to assign `ext.args` values to a specific process.

Both directives work within their defined scope: either for all tests within the `main.nf.test` file when written in the main `nextflow_process`, `nextflow_workflow` or `nextflow_pipeline` scope, or for a single test when written within the `test` scope.

### File paths

nf-test replaces paths in snapshots with a unique fingerprint (MD5 sum by default) to ensure file content consistency.

### Asserting channel item presence with `contains`

Use Groovy's `contains` and `collect` methods to assert the presence of items in channel output.

```groovy
// Example channel with tuples
def exampleChannel = [
  ['Bonjour', '/.nf-test/tests/c563c/work/65/b62f/Bonjour.json'],
  ['Hello', '/.nf-test/tests/c563c/work/65/fa20/Hello.json'],
  ['Hola', '/.nf-test/tests/c563c/work/65/85d0/Hola.json']
]

// Asserting a tuple's presence
testData = exampleChannel.collect { greeting, jsonPath -> [greeting, path(jsonPath).json] }
assert testData.contains(['Hello', path('./myTestData/Hello.json').json])

// Asserting a subset (greeting only)
testData = exampleChannel.collect { greeting, _ -> greeting }
assert testData.contains('Hello')
```

### Indexing

You can access elements in output channels using index notation, for example:

```groovy
log.get(0).get(1)
```

which is equivalent to

```groovy
log[0][1]
```

The first `get(0)` or `[0]` corresponds to the emitted channel object.

The `get(1)` or `[1]` corresponds to the second element of the channel object.
Most nf-core modules and pipelines emit two components: a meta map and the files or directories.
Specifying `get(1)` or `[1]` corresponds to the files or directories for recording in a snapshot.

### Debugging

You can use `println` statements to print variables during tests for debugging purposes.

Place print statements within the `then` block, before `assertAll`.

```nextflow
then {
    def unstable_patterns_auth = [
        '**/mapped_reads_gc-content_distribution.txt',
        '**/genome_gc_content_per_window.png',
        '**/*.{svg,pdf,html}',
        '*.{svg,pdf,html}',
        '**/DamageProfiler.log',
        ]

    println("unstable_patterns_auth: " + unstable_patterns_auth)

    assertAll(
        { assert snapshot( stable_content_authentication     , stable_name_authentication*.name   ).match("authentication") },
    ...
```

## Related documentation

The nf-test documentation provides additional information on:

- [Updating snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#updating-snapshots)
- [Cleaning obsolete snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#cleaning-obsolete-snapshots)
- [Constructing complex snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#constructing-complex-snapshots)

## Next steps

Learn how to write effective test assertions following nf-core guidelines:

- [Writing assertions](./assertions.md) covers nf-core guidelines and common assertion patterns
- [Advanced techniques](./advanced.md) covers complex patterns and troubleshooting
