---
title: "14. Error Testing & Edge Cases"
subtitle: Testing failure scenarios and edge cases
weight: 140
---

## Testing Expected Failures

### Invalid Input Testing

```groovy
test("Invalid file format") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/invalid_format.txt', checkIfExists: true)  // Not a FASTQ
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("Invalid file format") ||
                     process.stderr.contains("not in FASTQ format") }
        )
    }
}
```

### Parameter Validation Testing

```groovy
test("Invalid parameter values") {
    when {
        params {
            quality_threshold = -10  // Invalid: negative value
        }
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('test.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("Quality threshold must be positive") }
        )
    }
}

test("Missing required parameters") {
    when {
        params {
            // Missing required 'reference_genome' parameter
        }
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('test.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("Reference genome parameter is required") }
        )
    }
}
```

### Resource Limit Testing

```groovy
test("Memory limit exceeded") {
    when {
        params {
            max_memory = '100.MB'  // Very low memory limit
        }
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('large_dataset.bam', checkIfExists: true)  // Large file
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("OutOfMemoryError") ||
                     process.stderr.contains("exceeded memory limit") }
        )
    }
}
```

## Edge Case Testing

### Empty Input Files

```groovy
test("Empty input file") {
    setup {
        // Create empty file
        def emptyFile = new File("tests/empty.fastq")
        emptyFile.text = ""
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/empty.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            // Could succeed with empty output or fail gracefully
            { assert process.exitStatus in [0, 1] },
            {
                if (process.success) {
                    // Should produce empty or minimal output
                    assert process.out.size() >= 0
                } else {
                    // Should fail with appropriate message
                    assert process.stderr.contains("empty") ||
                           process.stderr.contains("no data")
                }
            }
        )
    }
}
```

### Single Record Input

```groovy
test("Single record input") {
    setup {
        // Create minimal FASTQ with single read
        def singleRead = new File("tests/single_read.fastq")
        singleRead.text = """@read1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
"""
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test', records: 1 ],
                file('tests/single_read.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match("single_record") }
        )
    }
}
```

### Maximum Size Input

```groovy
test("Maximum size input") {
    setup {
        // Create large test file (if not already exists)
        def largeFile = new File("tests/large_input.fastq")
        if (!largeFile.exists()) {
            largeFile.withWriter { writer ->
                (1..10000).each { i ->
                    writer.println("@read_${i}")
                    writer.println("A" * 150)
                    writer.println("+")
                    writer.println("I" * 150)
                }
            }
        }
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test', size: 'large' ],
                file('tests/large_input.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Verify reasonable execution time
            { assert process.trace.duration.toMinutes() < 30 },
            // Verify memory usage within limits
            { assert process.trace.peakRss < 8_000_000_000 }, // 8GB
            { assert snapshot(process.out).match("large_input") }
        )
    }
}
```

## Boundary Value Testing

### Numeric Boundaries

```groovy
[0, 1, 50, 99, 100].each { quality ->
    test("Quality threshold ${quality}") {
        when {
            params {
                quality_threshold = quality
            }
            process {
                """
                input[0] = [
                    [ id:'test', quality: ${quality} ],
                    file('test.fastq.gz', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match("quality_${quality}") }
            )
        }
    }
}
```

### String Length Boundaries

```groovy
['', 'a', 'a' * 255, 'a' * 256, 'a' * 1000].each { prefix ->
    test("Sample prefix length ${prefix.length()}") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test', prefix: '${prefix}' ],
                    file('test.fastq.gz', checkIfExists: true)
                ]
                """
            }
        }
        then {
            if (prefix.length() > 255) {
                assertAll(
                    { assert process.failed },
                    { assert process.stderr.contains("prefix too long") }
                )
            } else {
                assertAll(
                    { assert process.success },
                    { assert snapshot(process.out).match("prefix_${prefix.length()}") }
                )
            }
        }
    }
}
```

## Malformed Data Testing

### Corrupted File Headers

```groovy
test("Corrupted FASTQ header") {
    setup {
        def corruptedFastq = new File("tests/corrupted.fastq")
        corruptedFastq.text = """corrupted_header_no_at_symbol
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@read2
GCTAGCTAGCTAGCTA
+
JJJJJJJJJJJJJJJJ
"""
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/corrupted.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("invalid FASTQ format") ||
                     process.stderr.contains("malformed header") }
        )
    }
}
```

### Mismatched Sequence and Quality Lengths

```groovy
test("Mismatched sequence and quality lengths") {
    setup {
        def mismatchedFastq = new File("tests/mismatched.fastq")
        mismatchedFastq.text = """@read1
ATCGATCGATCGATCG
+
III
@read2
GCTA
+
JJJJJJJJJJJJJJJJ
"""
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/mismatched.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("sequence and quality length mismatch") ||
                     process.stderr.contains("length mismatch") }
        )
    }
}
```

## Concurrent Access Testing

### File Locking

```groovy
test("Concurrent file access") {
    when {
        process {
            """
            # Simulate concurrent access by multiple processes
            INPUT_FILE="shared_input.bam"

            # Try to acquire file lock
            exec 200>"\$INPUT_FILE.lock"
            if flock -n 200; then
                echo "File lock acquired"
                # Process the file
                samtools view "\$INPUT_FILE" | head -1000 > output1.sam
                # Release lock automatically when process exits
            else
                echo "File locked by another process, waiting..."
                flock 200  # Wait for lock
                samtools view "\$INPUT_FILE" | head -1000 > output1.sam
            fi

            input[0] = [
                [ id:'test' ],
                file('shared_input.bam', checkIfExists: true)
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

### Race Condition Testing

```groovy
test("Temporary file race conditions") {
    when {
        process {
            """
            # Create unique temporary files to avoid race conditions
            TEMP_PREFIX="temp_\${RANDOM}_\${task.index}"
            TEMP_FILE="\${TEMP_PREFIX}.tmp"

            # Ensure cleanup on exit
            trap "rm -f \${TEMP_PREFIX}*" EXIT

            # Use temporary file
            echo "Processing..." > "\$TEMP_FILE"
            sleep 1  # Simulate processing time
            mv "\$TEMP_FILE" "output.txt"

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

## Platform-Specific Testing

### Cross-Platform Path Handling

```groovy
test("Cross-platform path handling") {
    when {
        process {
            """
            # Handle paths correctly across platforms
            INPUT_PATH="${params.input_path}"

            # Normalize path separators
            if [[ "\$OSTYPE" == "msys" || "\$OSTYPE" == "cygwin" ]]; then
                # Windows/MSYS
                INPUT_PATH=\$(echo "\$INPUT_PATH" | sed 's|\\\\|/|g')
            fi

            # Verify path exists
            if [[ ! -f "\$INPUT_PATH" ]]; then
                echo "ERROR: Input file not found: \$INPUT_PATH"
                exit 1
            fi

            input[0] = [
                [ id:'test', platform: "\$OSTYPE" ],
                file("\$INPUT_PATH", checkIfExists: true)
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

### Line Ending Compatibility

```groovy
test("Line ending compatibility") {
    setup {
        // Create files with different line endings
        def unixFile = new File("tests/unix_endings.txt")
        unixFile.text = "line1\nline2\nline3\n"

        def windowsFile = new File("tests/windows_endings.txt")
        windowsFile.text = "line1\r\nline2\r\nline3\r\n"
    }

    when {
        process {
            """
            # Normalize line endings
            sed 's/\r\$//' tests/windows_endings.txt > normalized.txt

            input[0] = [
                [ id:'test' ],
                file('normalized.txt', checkIfExists: true)
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

## Timeout and Resource Exhaustion

### Process Timeout Testing

```groovy
test("Process timeout handling") {
    when {
        process {
            """
            # Set timeout for long-running process
            timeout 30s long_running_command input.txt output.txt || {
                echo "Process timed out, using fallback"
                cp input.txt output.txt  # Fallback action
                exit 0
            }

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
            { assert process.trace.duration.toSeconds() < 35 }, // Should complete within timeout + buffer
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Disk Space Testing

```groovy
test("Insufficient disk space") {
    when {
        process {
            """
            # Check available disk space
            AVAILABLE_SPACE=\$(df . | tail -1 | awk '{print \$4}')
            REQUIRED_SPACE=1000000  # 1GB in KB

            if [[ \$AVAILABLE_SPACE -lt \$REQUIRED_SPACE ]]; then
                echo "ERROR: Insufficient disk space"
                echo "Available: \$AVAILABLE_SPACE KB, Required: \$REQUIRED_SPACE KB"
                exit 1
            fi

            input[0] = [
                [ id:'test' ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        // Test should either succeed with sufficient space or fail gracefully
        assertAll(
            { assert process.exitStatus in [0, 1] },
            {
                if (process.failed) {
                    assert process.stderr.contains("Insufficient disk space")
                } else {
                    assert snapshot(process.out).match()
                }
            }
        )
    }
}
```

## Error Recovery Testing

### Graceful Degradation

```groovy
test("Graceful degradation") {
    when {
        process {
            """
            # Try optimal method, fall back to alternative
            if command -v fast_tool >/dev/null 2>&1; then
                fast_tool input.txt > output.txt 2>error.log || {
                    echo "Fast tool failed, trying alternative"
                    slow_but_reliable_tool input.txt > output.txt
                }
            else
                echo "Fast tool not available, using alternative"
                slow_but_reliable_tool input.txt > output.txt
            fi

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

### Retry Logic Testing

```groovy
test("Retry logic") {
    when {
        process {
            """
            # Retry logic for unreliable operations
            MAX_ATTEMPTS=3
            ATTEMPT=1

            while [[ \$ATTEMPT -le \$MAX_ATTEMPTS ]]; do
                echo "Attempt \$ATTEMPT of \$MAX_ATTEMPTS"

                if unreliable_command input.txt output.txt; then
                    echo "Command succeeded on attempt \$ATTEMPT"
                    break
                elif [[ \$ATTEMPT -eq \$MAX_ATTEMPTS ]]; then
                    echo "Command failed after \$MAX_ATTEMPTS attempts"
                    exit 1
                else
                    echo "Command failed, retrying in 5 seconds..."
                    sleep 5
                    ATTEMPT=\$((ATTEMPT + 1))
                fi
            done

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

## Test Data Validation

### Comprehensive Input Validation

```groovy
test("Input validation") {
    when {
        process {
            """
            INPUT_FILE="input.fastq.gz"

            # Validate file exists and is readable
            if [[ ! -r "\$INPUT_FILE" ]]; then
                echo "ERROR: Cannot read input file: \$INPUT_FILE"
                exit 1
            fi

            # Validate file is not empty
            if [[ ! -s "\$INPUT_FILE" ]]; then
                echo "ERROR: Input file is empty: \$INPUT_FILE"
                exit 1
            fi

            # Validate file format
            if ! zcat "\$INPUT_FILE" | head -4 | awk 'NR==1 {if(!\$0 ~ /^@/) exit 1} NR==3 {if(!\$0 ~ /^\+/) exit 1}'; then
                echo "ERROR: Invalid FASTQ format: \$INPUT_FILE"
                exit 1
            fi

            echo "Input validation passed"

            input[0] = [
                [ id:'test' ],
                file('input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.stdout.contains("Input validation passed") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Next Steps

Continue to [CI/CD Integration](./13_cicd_integration.md) to learn about integrating tests with continuous integration.
