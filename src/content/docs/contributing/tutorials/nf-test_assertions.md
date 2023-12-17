---
title: nf-test Assertions
description: A guide to using nf-test assertions for testing nf-core pipelines.
---

This document details various assertions used in nf-test for testing Nextflow pipelines. It serves as a guide for implementing effective testing strategies in pipeline development. For more information on nf-test, see the [nf-test documentation](https://code.askimed.com/nf-test/docs/getting-started/).

# Snapshots

Snapshots are used to compare the current output of a process, workflow, or function against a reference snapshot file (`*.nf.test.snap`).

## Using Snapshots

Create snapshots using the `snapshot` keyword. The `match` method checks if the snapshot corresponds to the expected data in the snap file. For example:

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

The first test run generates a json snapshot file. Subsequent runs compare against this file. Commit snapshot files with code changes and review them in your code review process.

## File Paths

nf-test replaces paths in snapshots with a unique fingerprint (md5 sum by default) to ensure file content consistency.

## Asserting the Presence of an Item in the Channel using `contains`

Groovy's `contains` and `collect` methods assert the presence of items in channel output.

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

## Indexing

You can access elements in output channels using index notation, for example:

```groovy
log.get(0).get(1)
```

which is equivalent to

```groovy
log[0][1]
```

## Additional Reading

- [Updating Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#updating-snapshots)
- [Cleaning Obsolete Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#cleaning-obsolete-snapshots)
- [Constructing Complex Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#constructing-complex-snapshots)

# nf-core guidelines for assertions

1. **Encapsulate Assertions in `assertAll()`**: Group all assertions within `assertAll()` for comprehensive testing.
2. **Minimum Requirement - Process Success + version.yml file**: Always check if the process completes successfully and make at least a snapshot of the version.yml.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out.version).match("version") }
)
```

3. **Capture as much as possible**: Best case scenario: make a snapshots to verify the complete output of your process. The absolute minimum is to check that the [output file exists](#file-exists-check), but try to check also for substrings, number of lines or similar.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

:::note
`process.out` will capture all output channels, both named and index based ones.
:::

## Additional cases:

4. **Handling Inconsistent md5sum**: Use specific content checks for elements with inconsistent md5sums.

5. **Module/Process Truth Verification**: Ensure snapshots accurately reflect the module/process functionality.

# Different Types of Assertions

## Simple & Straight-Forward

### Snapshot Entire Output Channel

_Motivation:_ Make sure all outputs are stable over changes.

```groovy {3}
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

_Explanation:_ Verifies process completion and output against a snapshot.

## Complex - Handling Inconsistent md5sum in Output Elements

### Snapshot a Specific Element in Output Channel

_Motivation:_ Create the snapshot for one specific output.

```groovy
assert snapshot(process.out.versions).match("versions")
```

_Explanation:_ Checks a specific element, in this case `versions`, in the output channel of a process against a predefined snapshot named "versions".

### File Exists Check

_Motivation:_ Snapshots of an output are unstable, i.e. they change between test runs, for example because they include a timestamp/file-path in the content.

- [BCLCONVERT](https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/tests/main.nf.test)

```groovy!
assert file(process.out.interop[0][1].find { file(it).name == "IndexMetricsOut.bin" }).exists()
```

_Explanation:_ Verifies the existence of a specific file, `IndexMetricsOut.bin`, in the output of a process.

### Snapshot Sorted List & Exclude a Specific File

_Motivation:_ I want to create a snapshot of different outputs, including several log files. I can't snapshot the whole output, because one file is changing between test runs.

- [BCLCONVERT](https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/tests/main.nf.test)

```groovy!
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

_Explanation:_ This creates a snapshot for all output files and of a sorted list from a log directory while excluding a specific file, `IndexMetricsOut.bin`, in the comparison. The existence of this excluded file is checked in the end.

### File Contains Check

```groovy {4} {title="bismark/align/tests/main.nf.test"}
with(process.out.report) {
    with(get(0)) {
        assert get(1).endsWith("hisat2_SE_report.txt")
        assert path(get(1)).readLines().last().contains("Bismark completed in")
    }
}
```

_Explanation:_ This checks if the last line of a report file contains a specific string and if the file name ends with "hisat2_SE_report.txt".

### Snapshot Selective Portion of a File

_Motivation:_ We can't make a snapshot of the whole file, because they are not stable, but we know a portion of the content should be stable, e.g. the timestamp is added in the 6th line, so we want to only snapshot the content of the first 5 lines.

```groovy
assert snapshot(file(process.out.aligned[0][1]).readLines()[0..4]).match()
```

_Explanation:_ Creates a snapshot of a specific portion (first five lines) of a file for comparison.

### Snapshot Selective Portion of a File & number of lines

_Motivation:_ We can't make a snapshot of the whole file, because they are not stable, but we know a portion of the content should be stable and the number of lines in it as well.

```groovy
def lines = path(process.out.file_out[0][1]).linesGzip
assertAll(
    { assert process.success },
    { assert snapshot(lines[0..5]).match("test_cat_zipped_zipped_lines") },
    { assert snapshot(lines.size()).match("test_cat_zipped_zipped_size") }
)
```

_Explanation:_ Verifies the content of the first six lines of a gzipped file, and the total number of lines in the file.

### ReadLines & Contains

_Motivation:_ We can't make a snapshot of the complete file, but we want to make sure that a specific substring is always present.

- [sratoolsncbisettings](https://github.com/nf-core/modules/blob/master/modules/nf-core/custom/sratoolsncbisettings/tests/main.nf.test)

```groovy
with(process.out.ncbi_settings) {
    assert path(get(0)).readLines().any { it.contains('/LIBS/GUID') }
    assert path(get(0)).readLines().any { it.contains('/libs/cloud/report_instance_identity') }
}
```

_Explanation:_ Checks if specific strings, `/LIBS/GUID` and `/libs/cloud/report_instance_identity` exist within the lines of an output file.

### Snapshot an Element in Tuple Output

_Motivation:_ We can't snapshot the whole tuple, but on element of the tuple has stable snapshots.

```groovy
assert snapshot(file(process.out.deletions[0][1])).match("deletions")
```

_Explanation:_ Validates an element within a tuple output against a snapshot.

### Snapshot Published File in Outdir

_Motivation:_ I want to check a specific file in the output is saved correctly and is stable between tests.

```groovy
params {
    outdir = "$outputDir"
}
```

```groovy
assert snapshot(path("$outputDir/kallisto/test/abundance.tsv")).match("abundance_tsv_single")
```

_Explanation:_ Confirms that a file saved in the specified output directory matches the expected snapshot.

### Assert File Name and Type

_Motivation:_ I don't know the exact location, know that at least the file type is fixed.

```groovy
assert process.out.classified_reads_fastq[0][1][0] ==~ ".*/test.classified_1.fastq.gz"
```

_Explanation:_ Ensures that a file from the output matches a specific pattern, indicating its type and name.

### Snapshot Selective File Names & Content

_Motivation:_ I want to include in the snapshot:

- the names of the files in `npa` & `npc` output channels
- The first line of the file in `npo` out channel
- The md5sum of the file in `npl` out channel

```groovy
assert snapshot(
    file(process.out.npa[0][1]).name,
    file(process.out.npc[0][1]).name,
    path(process.out.npo[0][1]).readLines()[0],
    path(process.out.npl[0][1])
).match()
```

_Explanation:_ Compares specific filenames and content of multiple files in a process output against predefined snapshots.

### Snapshot the Last 4 Lines of a Gzipped File in the gzip output channel

```groovy
path(process.out.gzip[0][1]).linesGzip[-4..-1]
```

_Explanation:_ Retrieves and allows the inspection of the last four lines of a gzipped file from the output channel.

### Assert a contains check in a gzipped file

_Motivation:_ I want to check the presence of a specific string or data pattern within a gzipped file

```groovy!
{ assert path(process.out.vcf[0][1]).linesGzip.toString().contains("MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242") }
```

_Explanation:_ check if a specific string (`"MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242"`) is present in the content of a gzipped file, specified by `path(process.out.vcf[0][1]).linesGzip.toString()`.

## Useful nf-test operators and functions

### Regular Expressions

The operator `==~` can be used to check if a string matches a regular expression:

```groovy
assert "/my/full/path/to/process/dir/example.vcf.pgen" ==~ ".*/example.vcf.pgen"
```

### Using `with()`

Instead of writing:

```groovy
assert process.out.imputed_plink2.size() == 1
        assert process.out.imputed_plink2[0][0] == "example.vcf"
        assert process.out.imputed_plink2[0][1] ==~ ".*/example.vcf.pgen"
        assert process.out.imputed_plink2[0][2] ==~ ".*/example.vcf.psam"
        assert process.out.imputed_plink2[0][3] ==~ ".*/example.vcf.pvar"
}
```

You can reduce redundancy using the `with()` command:

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

## Known Issues

When using nf-test in conjunction with container technologies like Docker, Singularity, or Conda, it's crucial to be aware of environment-specific issues that can arise, particularly regarding mismatched hashes. Here are some tips to handle such scenarios effectively:

### Tips for Handling Mismatched Hashes in Docker/Singularity/Conda

1. **Check for Consistent Environment Across Containers:**

   Ensure that the environment inside your Docker, Singularity, or Conda containers is consistent. Differences in installed packages, software versions, or underlying operating systems can lead to mismatched hashes.

2. **Use Identical Base Images:**

   When building Docker or Singularity containers, start from the same base image to minimize environmental differences. This consistency helps ensure that the software behaves the same across different executions.

3. **Pin Software Versions:**

   In your container definitions (Dockerfile, Singularity recipe, Conda environment file), explicitly pin software versions, including dependencies. This step reduces the chances of discrepancies due to updates or changes in the software.

4. **Isolate Non-Deterministic Elements:**

   Identify elements in your workflow that are inherently non-deterministic (such as timestamps or random number generation) and isolate them. Consider mocking these elements or designing your tests to accommodate such variability.

5. **Reproducibility in Conda Environments:**
   For Conda environments, use `conda list --explicit` to generate a list of all packages with their exact versions and builds. This approach ensures that you can recreate the identical environment later.

6. **Review Container Caching Mechanisms:**

   Be cautious with container caching mechanisms. Sometimes, cached layers in Docker might lead to using outdated versions of software or dependencies. Ensure that your caching strategy does not inadvertently introduce inconsistencies.

7. **Consistent Filesystem Paths:**

   Ensure that paths within the container and in the testing environment are consistent. Variations in paths can sometimes lead to unexpected behavior and hash mismatches.

8. **Regularly Update and Test:**
   Regularly update your containers and environment specifications, and re-run tests to ensure that everything continues to work as expected. This practice helps identify and resolve issues arising from environmental changes over time.

By following these tips, you can mitigate the risks of encountering mismatched hashes due to environment-specific issues in Docker, Singularity, and Conda when using nf-test for your Nextflow pipelines.
