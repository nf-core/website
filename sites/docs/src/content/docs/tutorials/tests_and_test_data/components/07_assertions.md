---
title: "7. nf-test Assertions"
subtitle: Comprehensive guide to nf-test assertions and verification patterns
weight: 70
---

This component covers various assertion patterns and techniques for effective testing with nf-test. Mastering these patterns is essential for creating robust and maintainable tests for nf-core pipelines and components.

## Snapshots

Snapshots are used to compare the current output of a process, workflow, or function against a reference snapshot file (`*.nf.test.snap`).

### Basic Snapshot Usage

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

The first test run generates a JSON snapshot file. Subsequent runs compare against this file. Commit snapshot files with code changes and review them in your code review process.

### Configuring Tests

nf-test allows specifying params or including config files:

```groovy
params {
    foo = 'bar'
}
```

```groovy
config "./nextflow.config"
```

Use withName selectors to assign `ext.args` values to a specific process. Both directives work within the scope they are defined in.

### File Path Handling

nf-test replaces paths in snapshots with a unique fingerprint (md5 sum by default) to ensure file content consistency.

## nf-core Guidelines for Assertions

### Core Requirements

1. **Encapsulate Assertions in `assertAll()`**: Group all assertions within `assertAll()` for comprehensive testing.

2. **Minimum Requirement**: Always check process success and version.yml file:

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out.versions).match("versions") }
)
```

3. **Capture as much as possible**: Best practice is to snapshot complete output:

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

:::note
`process.out` captures all output channels, both named and index-based ones.
:::

## Essential Assertion Patterns

### Simple & Straightforward

#### Snapshot Entire Output Channel

**Motivation:** Ensure all outputs are stable over changes.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

#### Snapshot Specific Element

**Motivation:** Create snapshot for one specific output.

```groovy
assert snapshot(process.out.versions).match("versions")
```

### File Verification Patterns

#### File Exists Check

**Motivation:** Snapshots are unstable due to timestamps/file-paths in content.

```groovy
assert file(process.out.interop[0][1].find { file(it).name == "IndexMetricsOut.bin" }).exists()
```

#### File Contains Check

**Motivation:** Verify specific content within files.

```groovy
with(process.out.report) {
    with(get(0)) {
        assert get(1).endsWith("hisat2_SE_report.txt")
        assert path(get(1)).readLines().last().contains("Bismark completed in")
    }
}
```

#### ReadLines & Contains

**Motivation:** Ensure specific substrings are always present.

```groovy
with(process.out.ncbi_settings) {
    assert path(get(0)).readLines().any { it.contains('/LIBS/GUID') }
    assert path(get(0)).readLines().any { it.contains('/libs/cloud/report_instance_identity') }
}
```

### Advanced Content Verification

#### Snapshot Selective File Portions

**Motivation:** Part of file content is stable, part is not (e.g., timestamps).

```groovy
// First 5 lines only
assert snapshot(file(process.out.aligned[0][1]).readLines()[0..4]).match()

// Lines and size of gzipped file
def lines = path(process.out.file_out[0][1]).linesGzip
assertAll(
    { assert process.success },
    { assert snapshot(lines[0..5]).match("test_cat_zipped_zipped_lines") },
    { assert snapshot(lines.size()).match("test_cat_zipped_zipped_size") }
)

// Last 4 lines of gzipped file
path(process.out.gzip[0][1]).linesGzip[-4..-1]
```

#### Assert Contains in Gzipped Files

**Motivation:** Check presence of specific data patterns in compressed files.

```groovy
{ assert path(process.out.vcf[0][1]).linesGzip.toString().contains("MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242") }
```

### Complex Output Handling

#### Snapshot Sorted List & Exclude Specific Files

**Motivation:** Snapshot multiple outputs while excluding unstable files.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(
        process.out.reports,
        process.out.versions,
        process.out.fastq,
        process.out.undetermined,
        file(process.out.logs.get(0).get(1)).list().sort(),
        process.out.interop.get(0).get(1).findAll { file(it).name != "IndexMetricsOut.bin" },
        ).match()
    },
    { assert file(process.out.interop.get(0).get(1).find { file(it).name == "IndexMetricsOut.bin" }).exists() }
)
```

#### Snapshot Element in Tuple Output

**Motivation:** Only one element of tuple has stable snapshots.

```groovy
assert snapshot(file(process.out.deletions[0][1])).match("deletions")
```

#### Snapshot Published File in Outdir

**Motivation:** Verify files saved correctly in output directory.

```groovy
params {
    outdir = "$outputDir"
}
```

```groovy
assert snapshot(path("$outputDir/kallisto/test/abundance.tsv")).match("abundance_tsv_single")
```

### Pattern Matching and Validation

#### Assert File Name and Type

**Motivation:** Verify file type when exact location is unknown.

```groovy
assert process.out.classified_reads_fastq[0][1][0] ==~ ".*/test.classified_1.fastq.gz"
```

#### Snapshot Selective File Names & Content

**Motivation:** Include specific aspects of multiple files in snapshot.

```groovy
assert snapshot(
    file(process.out.npa[0][1]).name,
    file(process.out.npc[0][1]).name,
    path(process.out.npo[0][1]).readLines()[0],
    path(process.out.npl[0][1])
).match()
```

### Handling Variable Directory Contents

#### Snapshotting Variable Files in Channel Directories

**Context:** Channel emits directory with mixed stable/unstable files.
**Motivation:** Snapshot stable files with md5sums, unstable files by name only.

```groovy
then {
    def stablefiles = []
    file(process.out.db.get(0).get(1)).eachFileRecurse{ file -> 
        if (!file.isDirectory() && !["database.log", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"].find {file.toString().endsWith(it)}) {
            stablefiles.add(file)
        } 
    }
    
    def unstablefiles = []
    file(process.out.db.get(0).get(1)).eachFileRecurse{ file -> 
        if (["database.log", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"].find {file.toString().endsWith(it)}) {
            unstablefiles.add(file.getName().toString())
        } 
    }

    assertAll(
        { assert process.success },
        { assert snapshot(
                stablefiles,
                unstablefiles,
                process.out.versions
            ).match() }
    )
}
```

:::note
Explicitly exclude directories in the first case because `eachFileRecurse` includes directories. If directories are included, nf-test looks inside and runs md5sum on files, potentially including files you wanted to exclude.
:::

### Binary File Verification

#### Snapshotting Variable Binary Files with File Size

**Context:** Binary files that cannot be asserted for valid contents.
**Motivation:** Check that binary file contains 'something' rather than just existence.

```groovy
// Exact file size (in bytes)
"malt/malt_index/ref.idx - correct file size: ${file("$outputDir/malt/malt_index/ref.idx").length()}",

// Minimum size (in bytes)
"malt/malt_index/ref.idx - minimum file size: ${file("$outputDir/malt/malt_index/ref.idx").length() >= 61616}",
```

## Channel Assertion Techniques

### Asserting Presence using `contains`

Use Groovy's `contains` and `collect` methods to assert presence of items in channel output:

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

Access elements in output channels using index notation:

```groovy
log.get(0).get(1)
// equivalent to
log[0][1]
```

- First `get(0)` or `[0]`: emitted channel object
- `get(1)` or `[1]`: second object (typically files/directories in nf-core modules)

## Useful nf-test Operators and Functions

### Regular Expressions

The operator `==~` checks if a string matches a regular expression:

```groovy
assert "/my/full/path/to/process/dir/example.vcf.pgen" ==~ ".*/example.vcf.pgen"
```

### Using `with()` for Cleaner Code

Instead of repetitive assertions:

```groovy
assert process.out.imputed_plink2.size() == 1
assert process.out.imputed_plink2[0][0] == "example.vcf"
assert process.out.imputed_plink2[0][1] ==~ ".*/example.vcf.pgen"
assert process.out.imputed_plink2[0][2] ==~ ".*/example.vcf.psam"
assert process.out.imputed_plink2[0][3] ==~ ".*/example.vcf.pvar"
```

Use `with()` to reduce redundancy:

```groovy
assert process.out.imputed_plink2
with(process.out.imputed_plink2) {
    assert size() == 1
    with(get(0)) {
        assert get(0) == "example.vcf"
        assert get(1) ==~ ".*/example.vcf.pgen"
        assert get(2) ==~ ".*/example.vcf.psam"
        assert get(3) ==~ ".*/example.vcf.pvar"
    }
}
```

## Debugging Techniques

### Using println for Debugging

Print statements can help debug test issues. They must go within the `then` block, prior to `assertAll`:

```groovy
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
        { assert snapshot( stable_content_authentication, stable_name_authentication*.name ).match("authentication") },
        // ... more assertions
    )
}
```

## Known Issues and Container Considerations

When using nf-test with Docker, Singularity, or Conda, be aware of environment-specific issues that can cause mismatched hashes:

### Tips for Handling Mismatched Hashes

1. **Consistent Environment**: Ensure consistent environments across containers
2. **Identical Base Images**: Use same base images for Docker/Singularity containers
3. **Pin Software Versions**: Explicitly pin software versions and dependencies
4. **Isolate Non-Deterministic Elements**: Identify and isolate non-deterministic elements
5. **Reproducible Conda Environments**: Use `conda list --explicit` for exact environment recreation
6. **Review Container Caching**: Be cautious with caching mechanisms
7. **Consistent Filesystem Paths**: Ensure path consistency within containers
8. **Regular Updates**: Regularly update and test containers

## Best Practices Summary

1. **Start Simple**: Begin with basic success checks and version snapshots
2. **Capture Everything Possible**: Aim to snapshot complete outputs when stable
3. **Handle Instability Gracefully**: Use selective snapshots for unstable content
4. **Use Meaningful Tags**: Name snapshots descriptively for easier debugging
5. **Debug Systematically**: Use println statements to understand output structure
6. **Review Snapshot Changes**: Always review snapshot file changes in code reviews
7. **Test Edge Cases**: Include tests for error conditions and boundary cases

## Next Steps

Continue to [Configuration Management](./07_configuration_management.md) to learn about managing test configurations and environments.

### Additional Reading

- [nf-test Documentation](https://code.askimed.com/nf-test/docs/getting-started/)
- [Updating Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#updating-snapshots)
- [Cleaning Obsolete Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#cleaning-obsolete-snapshots)
- [Constructing Complex Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#constructing-complex-snapshots) 