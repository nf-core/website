---
title: Advanced techniques
subtitle: Complex assertion patterns and troubleshooting
shortTitle: Advanced
---

Advanced assertion patterns handle complex testing scenarios where standard snapshots are insufficient.

## Snapshot sorted list and exclude specific files

**When to use:** You need to create snapshots of multiple outputs including log files, but one file changes between test runs.

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

This creates a snapshot for all output files and a sorted list from a log directory while excluding a specific file (`IndexMetricsOut.bin`) from the comparison. The existence of the excluded file is verified separately.

## Snapshot selective file portions

**When to use:** You cannot snapshot the whole file because it is unstable, but you know a portion of the content remains stable (for example, the timestamp is added after the first five lines).

```groovy
assert snapshot(file(process.out.aligned[0][1]).readLines()[0..4]).match()
```

This creates a snapshot of a specific portion (first five lines) of a file for comparison.

## Snapshot file portions with line count

**When to use:** You cannot snapshot the whole file because it is unstable, but you know a portion of the content and the total line count remain stable.

```groovy
def lines = path(process.out.file_out[0][1]).linesGzip
assertAll(
    { assert process.success },
    { assert snapshot(lines[0..5]).match("test_cat_zipped_zipped_lines") },
    { assert snapshot(lines.size()).match("test_cat_zipped_zipped_size") }
)
```

This verifies the content of the first six lines of a gzipped file and the total line count.

## Snapshot published files in output directory

**When to use:** You need to verify that a specific file in the output directory is saved correctly and remains stable between tests.

```groovy
params {
    outdir = "$outputDir"
}
```

```groovy
assert snapshot(path("$outputDir/kallisto/test/abundance.tsv")).match("abundance_tsv_single")
```

This confirms that a file saved in the specified output directory matches the expected snapshot.

## Snapshot selective filenames and content

**When to use:** You want to include in the snapshot:

- The filenames in `npa` and `npc` output channels
- The first line of the file in `npo` output channel
- The MD5 sum of the file in `npl` output channel

```groovy
assert snapshot(
    file(process.out.npa[0][1]).name,
    file(process.out.npc[0][1]).name,
    path(process.out.npo[0][1]).readLines()[0],
    path(process.out.npl[0][1])
).match()
```

This compares specific filenames and content of multiple files in a process output against predefined snapshots.

## Snapshot last lines of gzipped files

```groovy
path(process.out.gzip[0][1]).linesGzip[-4..-1]
```

This retrieves and allows inspection of the last four lines of a gzipped file from the output channel.

## Assert content in gzipped files

**When to use:** You need to verify the presence of a specific string or data pattern within a gzipped file.

```groovy!
{ assert path(process.out.vcf[0][1]).linesGzip.toString().contains("MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242") }
```

This checks if a specific string (`"MT192765.1\t10214\t.\tATTTAC\tATTAC\t29.8242"`) is present in the content of a gzipped file.

## Snapshot variable files in directory channels

**When to use:** A channel emits a directory and nf-test recursively lists all files, generating MD5 sums. Some files in the directory have unstable or empty MD5 sums. You need to snapshot all stable files with MD5 sums but only capture names of unstable files.

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

This code creates two lists of file paths within the emitted directory, filters them for stable and unstable files, and snapshots the lists.

The process:

1. Generate an empty list (`stablefiles`)
2. Retrieve the directory from the channel `db` using `get(1)` (rather than the meta)
3. Use `eachFileRecurse` to iterate through all files and directories
4. Append to the list (`.add(file)`) only files that are not directories and do not match the unstable filenames (`"database.log"`, `"database.fastaid2LCAtaxid"`, `"database.taxids_with_multiple_offspring"`)
5. Create a second list (`unstablefiles`) and append only files that match the unstable filenames
6. For unstable files, append only the filename (not the full path) by converting to a string (`getName().toString()`)

Capture both lists in the snapshot using `assert snapshot().match()`.

:::note
Explicitly exclude directories because `eachFileRecurse` includes them when listing files.

If directories are included in the list of files to snapshot, nf-test looks inside any directory in the list (such as `stablefiles`) and runs an MD5 sum on any file in that directory.
Even if you explicitly exclude the file during the `eachFileRecurse` and `find` functions, nf-test will still process it via the directory.

Excluding directories prevents accidental duplicate listing of files you want to exclude.
:::

## Snapshot variable binary files with file size

**When to use:** Your tool produces a binary file that cannot be validated by checking string contents, but you want to verify the binary file contains valid content rather than checking only for file existence.

Compare an exact file size (in bytes):

```nextflow
"malt/malt_index/ref.idx - correct file size: ${file("$outputDir/malt/malt_index/ref.idx").length()}",
```

Check for a minimum size (in bytes):

```nextflow
"malt/malt_index/ref.idx - minimum file size: ${file("$outputDir/malt/malt_index/ref.idx").length() >= 61616}",
```

Binary files with variable contents cannot use MD5 sums because the hash differs each time. As binary files, you cannot search for plain text strings to verify specific content.

Checking only for file existence may be insufficient, as some tools can produce binary files with inadequate contents. If your tool produces a binary file with stable size (despite content variability), or if a working binary file exceeds a particular size, you can use the file size (in bytes) to assert validity.

## Useful operators and functions

### Regular expressions

Use the `==~` operator to check if a string matches a regular expression:

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

## Troubleshooting

When using nf-test with container technologies like Docker, Singularity, or Conda, environment-specific issues can arise, particularly mismatched hashes.

How tp handle mismatched hashes in containers:

1. **Check for consistent environments across containers**

   Ensure that the environment inside your Docker, Singularity, or Conda containers is consistent. Differences in installed packages, software versions, or underlying operating systems can lead to mismatched hashes.

2. **Use identical base images**

   When building Docker or Singularity containers, start from the same base image to minimise environmental differences. This consistency helps ensure that the software behaves the same across different executions.

3. **Pin software versions**

   In your container definitions (Dockerfile, Singularity recipe, Conda environment file), explicitly pin software versions, including dependencies. This reduces the chances of discrepancies due to updates or changes in the software.

4. **Isolate non-deterministic elements**

   Identify elements in your workflow that are inherently non-deterministic (such as timestamps or random number generation) and isolate them. Consider mocking these elements or designing your tests to accommodate such variability.

5. **Ensure reproducibility in Conda environments**

   For Conda environments, use `conda list --explicit` to generate a list of all packages with their exact versions and builds. This ensures that you can recreate the identical environment later.

6. **Review container caching mechanisms**

   Be cautious with container caching mechanisms. Cached layers in Docker might lead to using outdated versions of software or dependencies. Ensure that your caching strategy does not inadvertently introduce inconsistencies.

7. **Maintain consistent file system paths**

   Ensure that paths within the container and in the testing environment are consistent. Variations in paths can lead to unexpected behaviour and hash mismatches.

8. **Regularly update and test**

   Regularly update your containers and environment specifications, and re-run tests to ensure that everything continues to work as expected. This practice helps identify and resolve issues arising from environmental changes over time.
