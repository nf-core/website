---
title: "9. Snapshot Management"
subtitle: Working with snapshots and test outputs
weight: 90
---

## Understanding Snapshots

Snapshots capture the exact output of your tests and store them for comparison in future test runs. They provide a reliable way to detect changes in module, subworkflow, or pipeline outputs.

## Basic Snapshot Usage

### Creating Snapshots

```groovy
test("Basic snapshot test") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Updating Snapshots

```bash
# Update all snapshots
nf-test test --update-snapshot

# Update specific test snapshots
nf-test test modules/nf-core/fastqc/tests/main.nf.test --update-snapshot
```

## Snapshot Structure

### Snapshot Files

Snapshots are stored in `.nf-test/snapshots/` directory:

```
.nf-test/snapshots/
├── modules/
│   └── nf-core/
│       └── fastqc/
│           └── tests/
│               └── main.nf.test.snap
└── subworkflows/
    └── nf-core/
        └── bam_sort_stats_samtools/
            └── tests/
                └── main.nf.test.snap
```

### Snapshot Content

```json
{
  "Single-end reads": {
    "content": [
      {
        "0": [
          [
            {
              "id": "test",
              "single_end": true
            },
            "test.html:md5sum:123abc456def789"
          ]
        ],
        "1": [
          [
            {
              "id": "test",
              "single_end": true
            },
            "test_fastqc.zip:md5sum:987fed654cba321"
          ]
        ],
        "2": [
          "versions.yml:md5sum:abc123def456"
        ]
      }
    ],
    "timestamp": "2024-01-15T10:30:45.123456"
  }
}
```

## Advanced Snapshot Patterns

### Selective Snapshots

```groovy
then {
    assertAll(
        { assert process.success },
        // Snapshot only specific outputs
        { assert snapshot(
            process.out.html,
            process.out.zip
        ).match() }
    )
}
```

### Named Snapshots

```groovy
test("Multiple conditions") {
    when {
        params {
            condition = 'stringent'
        }
        process {
            """
            input[0] = [
                [ id:'test', condition: params.condition ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match("stringent_condition") }
        )
    }
}
```

### Conditional Snapshots

```groovy
then {
    assertAll(
        { assert process.success },
        {
            if (params.generate_report) {
                assert snapshot(process.out.report).match("with_report")
            } else {
                assert snapshot(process.out).match("no_report")
            }
        }
    )
}
```

## Snapshot Filtering

### Excluding Variable Content

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Filter out timestamps and other variable content
            process.out.log.collect { 
                it.readLines().findAll { line -> 
                    !line.contains("timestamp") && !line.contains("execution_id")
                }
            }
        ).match() }
    )
}
```

### Content Transformation

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Transform output before snapshotting
            process.out.stats.collect { meta, stats ->
                [
                    meta: [id: meta.id], // Keep only stable meta fields
                    stats: stats.readLines()[0..10] // Only first 10 lines
                ]
            }
        ).match() }
    )
}
```

## File Content Snapshots

### Text File Content

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Snapshot file content
            process.out.txt.collect { meta, txt -> txt.text }
        ).match() }
    )
}
```

### Structured Data

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Parse and snapshot JSON content
            process.out.json.collect { meta, json ->
                new groovy.json.JsonSlurper().parse(json)
            }
        ).match() }
    )
}
```

## Binary File Handling

### Checksum Snapshots

```groovy
then {
    assertAll(
        { assert process.success },
        // Binary files automatically use checksums
        { assert snapshot(process.out.bam).match() },
        // Explicit checksum for verification
        { assert snapshot(
            process.out.bam.collect { meta, bam ->
                [meta: meta, checksum: bam.md5()]
            }
        ).match("bam_checksums") }
    )
}
```

### Size Verification

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            process.out.bam.collect { meta, bam ->
                [
                    id: meta.id,
                    size: bam.size(),
                    exists: bam.exists()
                ]
            }
        ).match("bam_properties") }
    )
}
```

## Snapshot Best Practices

### 1. Meaningful Test Names

```groovy
test("Single-end reads with quality filtering") {
    when {
        params {
            min_quality = 30
        }
        // ... test implementation
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match("single_end_quality_30") }
        )
    }
}
```

### 2. Stable Content Only

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Remove timestamps and process IDs
            process.out.log.collect { meta, log ->
                log.readLines()
                   .findAll { !it.contains("timestamp") }
                   .findAll { !it.contains("process_id") }
                   .join('\n')
            }
        ).match() }
    )
}
```

### 3. Hierarchical Organization

```groovy
test("Complex workflow") {
    then {
        assertAll(
            { assert workflow.success },
            // Separate snapshots for different components
            { assert snapshot(workflow.out.qc_reports).match("qc_reports") },
            { assert snapshot(workflow.out.alignments).match("alignments") },
            { assert snapshot(workflow.out.variants).match("variants") }
        )
    }
}
```

## Snapshot Maintenance

### Updating Selective Snapshots

```bash
# Update snapshots for specific tag
nf-test test --tag fastqc --update-snapshot

# Update snapshots for specific pattern
nf-test test "modules/nf-core/*/tests/main.nf.test" --update-snapshot
```

### Viewing Snapshot Differences

```bash
# View what changed (requires git)
git diff .nf-test/snapshots/

# Compare snapshots before updating
nf-test test --verbose  # Shows differences
```

### Cleaning Old Snapshots

```bash
# Clean unused snapshots
nf-test clean

# Remove all snapshots (careful!)
rm -rf .nf-test/snapshots/
```

## Snapshot Validation

### Pre-commit Checks

```groovy
test("Snapshot validation") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Validate before snapshotting
            {
                process.out.txt.each { meta, txt ->
                    assert txt.exists()
                    assert txt.size() > 0
                    assert txt.text.contains("expected_content")
                }
            },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Cross-platform Compatibility

```groovy
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            // Normalize line endings for cross-platform compatibility
            process.out.txt.collect { meta, txt ->
                txt.text.replaceAll(/\r\n/, '\n')
            }
        ).match() }
    )
}
```

## Troubleshooting Snapshots

### Common Issues

1. **Snapshot Mismatches**: Often due to timestamps or system-specific paths
2. **Binary File Differences**: Check if files are truly different or just metadata
3. **Missing Snapshots**: Run with `--update-snapshot` for new tests

### Debug Snapshot Content

```groovy
then {
    assertAll(
        { assert process.success },
        {
            // Debug: print actual content before snapshotting
            process.out.each { channel ->
                println "Channel content: ${channel}"
            }
        },
        { assert snapshot(process.out).match() }
    )
}
```

## Next Steps

Continue to [Advanced Testing Patterns](./10_advanced_testing_patterns.md) to learn complex testing scenarios. 