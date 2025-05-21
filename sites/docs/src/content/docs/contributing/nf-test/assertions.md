---
title: "nf-test: Example assertions"
subtitle: A guide to using nf-test assertions for testing nf-core pipelines.
parentWeight: 20
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

## Assigning parameters or configs

nf-test allows to specify params or including config files.

```groovy
params {
  foo = 'bar'
}
```

```groovy
config "./nextflow.config"
```

Use withName selectors to assign `ext.args` values to a specific process.

Both these directives work within the scope they are defined in.
So either for the full set of test within the `main.nf.test` file if written in the main `nextflow_process`, `nextflow_workflow` or `nextflow_pipeline` scope, or for a single test if written within the `test` scope.

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
process.out.log.get(0).get(1)
```

which is equivalent to

```groovy
process.out.log[0][1]
```

or even better to be more descriptive using the collect operator as

```groovy
process.out.log.collect { meta, log -> log }
```

The first `get(0)` or `[0]` corresponds to the emitted channel object itself.

The `get(1)` or `[1]` corresponds to the second object of the channel object.
Most nf-core modules and pipelines typically emit two sub-components of an object: a meta map and the file(s)/directories etc.
Specifying `get(q)` or `[1]` thus corresponds to the file(s)/directories for recording in a snapshot.

## Debugging

When you assign variables that you inject into the `assertAll`, you can use `println` statements to print these variables during the test for debugging purposes.

The print statements must go within the `then` block, and prior `assertAll`.

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
    { assert snapshot(process.out.versions).match() }
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

4. **Configuration of `ext.args` in tests**: Module nf-tests SHOULD use a single nextflow.config to supply ext.args to a module. They can be defined in the when block of a test under the params scope.

```groovy {title="main.nf.test"}
config './nextflow.config'

when {
  params {
    module_args = '--extra_opt1 --extra_opt2'
  }
  process {
    """
    input[0] = [
      [ id:'test1', single_end:false ], // meta map
      file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)
    ]
    """
  }
}
```

```groovy {title="nextflow.config"}
process {
  withName: 'MODULE' {
    ext.args = params.module_args
  }
}
```

5. **Include a stub test and capture versions directly**: Each module should include a stub test. You can directly capture the content of the versions YAML file in the snapshot for better readability.

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(
        process.out,
        path(process.out.versions[0]).yaml
    ).match() }
)
```


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

_Explanation:_ Checks a specific element, in this case `versions`, in the output channel of a process against a predefined snapshot element named "versions".

### File Exists and Matches Name Check

_Motivation:_ Snapshots of an output are unstable, i.e. they change between test runs, for example because they include a timestamp/file-path in the content.

- [BCLCONVERT](https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/tests/main.nf.test)

```groovy
{ assert process.out.interop.collect { meta, interop ->
                                        interop.find { interopfile -> file(interopfile).name == "IndexMetricsOut.bin" } } }
```

_Explanation:_ Verifies the existence of a specific file, `IndexMetricsOut.bin`, in the output of a process.

### Snapshot Sorted List & Exclude a Specific File

_Motivation:_ I want to create a snapshot of different outputs, including several log files. I can't snapshot the whole output, because one file is changing between test runs.

- [BCLCONVERT](https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/tests/main.nf.test)

```groovy
assertAll(
    { assert process.success },
    { assert snapshot(
        process.out.fastq,
        process.out.fastq_idx,
        process.out.undetermined.collect { meta, fastq -> file(fastq).name },
        process.out.undetermined_idx,
        process.out.reports,
        process.out.logs.collect { meta, logs -> file(logs).list().sort() },
        process.out.interop.collect { meta, interop ->
                                    interop.findAll { interopfile ->
                                    file(interopfile).name != "IndexMetricsOut.bin" } },
        process.out.versions
    ).match() },
    { assert process.out.interop.collect { meta, interop ->
                                                        interop.find { interopfile -> file(interopfile).name == "IndexMetricsOut.bin" } } }
)
```

_Explanation:_ This creates a snapshot for all output files and of a sorted list from a log directory while excluding a specific file, `IndexMetricsOut.bin`, in the comparison. The name of this excluded file is instead included in the snapshot.

### File Contains Check

- [BISMARK/ALIGN](https://github.com/nf-core/modules/blob/master/modules/nf-core/bismark/align/tests/main.nf.test)

```groovy
process.out.report.collect { meta, report ->
                             file(report).readLines().contains("Number of alignments with a unique best hit from the different alignments:\t5009")
                           },
```

_Explanation:_ This checks if the last line of a report file contains a specific string.

### Snapshot Selective Portion of a File

_Motivation:_ We can't make a snapshot of the whole file, because they are not stable, but we know a portion of the content should be stable, e.g. the timestamp is added in the 6th line, so we want to only snapshot the content of the first 5 lines.

```groovy
assert snapshot(
          process.out.aligned.collect { meta, aligned ->
                                        file(aligned).readLines()[0..4]
                                      }
).match()
```

_Explanation:_ Creates a snapshot of a specific portion (first five lines) of a file for comparison.

### Snapshot Selective Portion of a File & number of lines

_Motivation:_ We can't make a snapshot of the whole file, because they are not stable, but we know a portion of the content should be stable and the number of lines in it as well.

```groovy
def lines = process.out.file_out.collect { meta, outfile -> path(outfile).linesGzip }
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
assert snapshot(process.out.deletions.collect { meta, deletions -> deletions }).match("deletions")
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
    process.out.npa.collect { meta, npa -> file(npa).name },
    process.out.npc.collect { meta, npc -> file(npc).name },
    process.out.npo.collect { meta, npo -> file(npo).readLines()[0] },
    process.out.npl.collect { meta, npl -> file(npl) }
).match()
```

_Explanation:_ Compares specific filenames and content of multiple files in a process output against predefined snapshots.

### Snapshot the Last 4 Lines of a Gzipped File in the gzip output channel

```groovy
assert snapshot(
    process.out.gzip.collect { meta, gzip -> gzip.linesGzip[-4..-1] }
)
```

_Explanation:_ Retrieves and allows the inspection of the last four lines of a gzipped file from the output channel.

### Assert a contains check in a gzipped file

_Motivation:_ I want to check the presence of a specific string or data pattern within a gzipped file

```groovy
assert process.out.vcf.collect { meta, vcf ->
                                  path(vcf).linesGzip.toString().contains("MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242")
}
```

_Explanation:_ check if a specific string (`"MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242"`) is present in the content of a gzipped file, specified by `path(process.out.vcf[0][1]).linesGzip.toString()`.

### Snapshotting variable files in a channel emitting a directory

_Context_: If a channel emits just a directory, by default nf-test will recursively list all files in that and all sub directories, and generate md5sums of all the files.
However, in some cases, _some_ of the files in the directory may have unstable/empty md5sums.
I want to snapshot all stable files with md5sums, but only snapshot names of unstable files.
_Motivation_: I want to snapshot all files with stable md5sums, but only snapshot names of unstable files.

```bash
then {
    def stablefiles = []
    file(process.out.db.get(0).get(1)).eachFileRecurse{ file -> if (!file.isDirectory() && !["database.log", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"].find {file.toString().endsWith(it)}) {stablefiles.add(file)} }
    def unstablefiles = []
    file(process.out.db.get(0).get(1)).eachFileRecurse{ file -> if (["database.log", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"].find {file.toString().endsWith(it)}) {unstablefiles.add(file.getName().toString())} }

    assertAll(
        { assert process.success },
        { assert snapshot(
                stableFiles,
                stableNames,
                process.out.versions
            ).match() }
    )
}
```

_Explanation_: We create two lists of files paths within the emitted directory, filter these two for stable and unstable files respectively, and snapshot the lists of paths.

In more detail, we generate an empty list (`stablefiles`). We then retrieve the directory from the channel `db` using `get(1)` (rather than the meta), and retrieve all files and directories that are inside that directory using `endFileRecurse` and, however we only append to the list (`.add(file)`) those files that are not a directory and not paths that end in (`endsWith`) the file names identified as unstable (`"database.log", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"`).

We then do the reverse (`unstablefiles`), where we loop again through the directory, but this time append only files that _do_ match the identified unstable file names.
However do not append the path itself, but just the filename by converting to a string (`getName().toString()`) when adding to the list.

These two lists of stable paths and unstable names can be captured in the snapshot in an `assert snapshot().match()`.

:::note
We have to explicitly exclude directories in the first case, because `eachFileRecurse` includes directories when listing all files.

If directories are included in the list of files to be snapshot, nf-test by default looks inside any directory in the list (here called `stablefiles`) and also runs an md5sum on any file in the listed directory.
Therefore, even if you explicitly exclude the file during the `endFileRercurse` and `find` function, and thus it is not explicitly in the `stablefiles` list itself, the file will still be picked by nf-test via the directory.

Therefore, by excluding directories, you do not get an accidental 'double' listing of files you wish to exclude.
:::

### Snapshotting variable binary files with file size

_Context_: You have a tool that always produces a binary file that cannot be asserted for valid contents such as a string.

_Motivation_: You want to be able to still check that the binary file contains 'something' rather than just the existence of the file.

To compare an exact file size (in bytes)

```nextflow
"malt/malt_index/ref.idx - correct file size: ${file("$outputDir/malt/malt_index/ref.idx").length()}",
```

To check for a minimum size (in bytes)

```nextflow
"malt/malt_index/ref.idx - minimum file size: ${file("$outputDir/malt/malt_index/ref.idx").length() >= 61616}",
```

_Explanation_: When you have a binary file that can have variable contents, you cannot use a md5sum, as the md5sum hash will be different each time. Then, as it is a binary file, you cannot easily search for plain text strings to check that the specific string is present in the file.

While you could check simply for the existence of a file, it may be that some tools can produce binary files that have 'insufficent' contents for it to work.
If you know that your tool produces a binary file _size_ that is stable (despite variability), or you know that a 'working' binary file exceeds a particular size, you can use the file size (in bytes) to assert the file is 'valid'.

## Useful nf-test operators and functions

### Regular Expressions

The operator `==~` can be used to check if a string matches a regular expression:

```groovy
assert "/my/full/path/to/process/dir/example.vcf.pgen" ==~ ".*/example.vcf.pgen"
```

## Known Issues

When using nf-test in conjunction with container technologies like Docker, Singularity, or Conda, it's crucial to be aware of environment-specific issues that can arise, particularly regarding mismatched hashes. Here are some tips to handle such scenarios effectively:

### Tips for Handling Mismatched Hashes in Docker/Singularity/Conda

1. **Check for Consistent Environment Across Containers:**

   Ensure that the environment inside your Docker, Singularity, or Conda containers is consistent. Differences in installed packages, software versions, or underlying operating systems can lead to mismatched hashes.

2. **Use Identical Base Images:**

   When building Docker or Singularity containers, start from the same base image to minimize environmental differences. This consistency helps ensure that the software behaves the same across different executions.

3. **Pin Software Versions:**

   In your container definitions (Dockerfile, Singularity recipe, Conda environment file), explicitly pin software versions, including dependencies. This practice reduces the likelihood of discrepancies due to unexpected updates or changes in software behavior.

4. **Isolate Non-Deterministic Elements:**

   Identify and isolate elements in your workflow that are inherently non-deterministic (such as timestamps, random number generation, or file paths). Consider mocking these elements, using fixed seeds for random processes, or designing your tests to accommodate such variability through pattern matching or partial comparisons.

5. **Ensure Reproducibility in Conda Environments:**
   For Conda environments, use `conda list --explicit > environment.lock.yml` to generate a detailed list of all packages with their exact versions and builds. This approach ensures that you can recreate the identical environment later, minimizing test inconsistencies.

6. **Review Container Caching Mechanisms:**

   Be cautious with container caching mechanisms. Cached layers in Docker might lead to using outdated versions of software or dependencies. Consider using the `--no-cache` flag when building containers for testing, or implement a consistent caching strategy that doesn't compromise reproducibility.

7. **Maintain Consistent Filesystem Paths:**

   Ensure that paths within the container and in the testing environment are consistent. Absolute vs. relative paths or variations in path structure can lead to unexpected behavior and hash mismatches. When possible, use path normalization in your test assertions.

8. **Implement Regular Updates and Testing:**
   Regularly update your containers and environment specifications, and re-run tests to ensure continued functionality. Schedule periodic reviews of your test environment to identify and resolve issues arising from environmental changes, dependency updates, or modifications to underlying infrastructure.

By following these guidelines, you can significantly reduce the risks of encountering mismatched hashes and other environment-specific issues when using nf-test with Docker, Singularity, and Conda in your Nextflow pipelines.
