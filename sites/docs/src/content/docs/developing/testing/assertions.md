---
title: nf-core guidelines for assertions
subtitle: nf-core guidelines and common patterns
shortTitle: Assertions
---

nf-core guidelines for assertions:

1. **Encapsulate assertions in `assertAll()`**: Group all assertions within `assertAll()` for comprehensive testing.
2. **Minimum requirement - process success and version.yml file**: Always check if the process completes successfully and make at least a snapshot of the `version.yml`.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out.versions).match("versions") }
)
```

3. **Capture as much as possible**: Make snapshots to verify the complete output of your process. At minimum, check that the [output file exists](#file-exists-check), but also verify substrings, line counts, or similar properties.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

:::note
`process.out` captures all output channels, both named and index-based.
:::

4. **Handle inconsistent MD5 sums**: Use specific content checks for elements with inconsistent MD5 sums.

5. **Verify module and process functionality**: Ensure snapshots accurately reflect the module or process functionality.

## Assertion patterns

### Simple patterns

#### Snapshot entire output channel

**When to use:** Verify that all outputs remain stable through changes.

```groovy {3}
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

This verifies process completion and output against a snapshot.

### Complex patterns

These patterns handle cases where outputs have inconsistent MD5 sums or require selective snapshotting.

#### Snapshot a specific output element

**When to use:** Create the snapshot for one specific output.

```groovy
assert snapshot(process.out.versions).match("versions")
```

This checks a specific element (in this case `versions`) in the output channel of a process against a predefined snapshot named `versions`.

#### File exists check

**When to use:** Snapshots of an output are unstable, changing between test runs because they include timestamps or file paths in the content.

- [BCLCONVERT](https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/tests/main.nf.test)

```groovy!
assert file(process.out.interop[0][1].find { file(it).name == "IndexMetricsOut.bin" }).exists()
```

This verifies the existence of a specific file (`IndexMetricsOut.bin`) in the process output.

#### File contains check

```groovy {4} {title="bismark/align/tests/main.nf.test"}
with(process.out.report) {
    with(get(0)) {
        assert get(1).endsWith("hisat2_SE_report.txt")
        assert path(get(1)).readLines().last().contains("Bismark completed in")
    }
}
```

This checks if the last line of a report file contains a specific string and if the filename ends with `hisat2_SE_report.txt`.

#### Check file content with readLines

**When to use:** You cannot snapshot the complete file but need to verify that specific substrings are always present.

- [sratoolsncbisettings](https://github.com/nf-core/modules/blob/master/modules/nf-core/custom/sratoolsncbisettings/tests/main.nf.test)

```groovy
with(process.out.ncbi_settings) {
    assert path(get(0)).readLines().any { it.contains('/LIBS/GUID') }
    assert path(get(0)).readLines().any { it.contains('/libs/cloud/report_instance_identity') }
}
```

This checks if specific strings (`/LIBS/GUID` and `/libs/cloud/report_instance_identity`) exist within the lines of an output file.

#### Snapshot tuple element

**When to use:** You cannot snapshot the whole tuple, but one element has stable snapshots.

```groovy
assert snapshot(file(process.out.deletions[0][1])).match("deletions")
```

This validates an element within a tuple output against a snapshot.

#### Assert filename and type

**When to use:** You don't know the exact file location but the file type is fixed.

```groovy
assert process.out.classified_reads_fastq[0][1][0] ==~ ".*/test.classified_1.fastq.gz"
```

This ensures that a file from the output matches a specific pattern, indicating its type and name.

## Next steps

[Advanced techniques](./advanced.md) covers complex assertion patterns for handling edge cases, useful operators, and troubleshooting guidance.
