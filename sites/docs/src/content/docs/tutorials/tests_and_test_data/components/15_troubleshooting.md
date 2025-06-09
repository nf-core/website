---
title: "17. Troubleshooting"
subtitle: Common issues and solutions
weight: 170
---

## Common Test Failures

### Snapshot Mismatches

**Problem**: Tests fail with snapshot differences

```
FAILED - snapshot does not match:
Expected: test.html:md5:abc123
Actual:   test.html:md5:def456
```

**Solutions**:

1. **Check for timestamps or variable content**:

   ```groovy
   // Filter variable content before snapshotting
   { assert snapshot(
       process.out.html.collect { meta, html ->
           html.text.replaceAll(/\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}/, 'TIMESTAMP')
       }
   ).match() }
   ```

2. **Update snapshots if changes are expected**:

   ```bash
   nf-test test --update-snapshot
   ```

3. **Check file permissions and paths**:
   ```bash
   # Verify files are accessible
   ls -la .nf-test/snapshots/
   ```

### File Not Found Errors

**Problem**: `checkIfExists: true` fails

```
ERROR: File does not exist: /path/to/test_data.fastq.gz
```

**Solutions**:

1. **Check test data paths**:

   ```groovy
   // Verify path construction
   def testFile = file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz')
   println "Looking for file at: ${testFile.toAbsolutePath()}"
   ```

2. **Validate test data parameter**:

   ```groovy
   // Add to .nf-test.config
   params {
       modules_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
   }
   ```

3. **Use local test data as fallback**:
   ```groovy
   def testFile = file('tests/data/local_test.fastq.gz').exists() ?
                  file('tests/data/local_test.fastq.gz') :
                  file(params.modules_testdata_base_path + 'path/to/remote_test.fastq.gz', checkIfExists: true)
   ```

### Container Issues

**Problem**: Container not found or access denied

```
ERROR: Docker image 'biocontainers/fastqc:0.11.9--0' not found
```

**Solutions**:

1. **Check container availability**:

   ```bash
   docker pull biocontainers/fastqc:0.11.9--0
   singularity pull docker://biocontainers/fastqc:0.11.9--0
   ```

2. **Use alternative registries**:

   ```groovy
   // In nextflow.config
   process {
       withName: 'FASTQC' {
           container = 'quay.io/biocontainers/fastqc:0.11.9--0'
       }
   }
   ```

3. **Test without containers first**:
   ```bash
   nf-test test --profile conda
   ```

### Memory/Resource Issues

**Problem**: Out of memory or resource limits exceeded

```
ERROR: Process 'TOOL' terminated with exit status 137 (SIGKILL - Out of memory)
```

**Solutions**:

1. **Increase resource limits**:

   ```groovy
   // In test nextflow.config
   process {
       withName: 'TOOL' {
           memory = '8.GB'
           cpus = 4
       }
   }
   ```

2. **Use smaller test data**:

   ```groovy
   // Create minimal test dataset
   setup {
       def smallFile = new File("tests/small_input.fastq")
       smallFile.text = """@read1
   ATCGATCGATCGATCG
   +
   IIIIIIIIIIIIIIII
   """
   }
   ```

3. **Monitor resource usage**:
   ```groovy
   then {
       assertAll(
           { assert process.success },
           { assert process.trace.peakRss < 2_000_000_000 }, // 2GB limit
           { assert snapshot(process.out).match() }
       )
   }
   ```

## Debugging Strategies

### Verbose Output

Enable detailed logging:

```bash
# Maximum verbosity
nf-test test --verbose --debug

# Nextflow-specific debugging
nf-test test -with-trace -with-report -with-timeline
```

### Work Directory Inspection

```bash
# Find work directories for failed processes
find work -name '.exitcode' -exec dirname {} \;

# Examine failed process
ls -la work/ab/cd123456/
cat work/ab/cd123456/.command.sh
cat work/ab/cd123456/.command.out
cat work/ab/cd123456/.command.err
```

### Interactive Debugging

```groovy
test("Debug test") {
    when {
        process {
            """
            # Add debugging output
            echo "DEBUG: Starting process with inputs:"
            echo "Input file: \$1"
            ls -la \$1

            # Check environment
            echo "DEBUG: Environment variables:"
            env | grep -E "(PATH|JAVA|NEXTFLOW)" | sort

            # Tool version check
            echo "DEBUG: Tool version:"
            my_tool --version

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
            // Print debug information
            { println "STDOUT: ${process.stdout}" },
            { println "STDERR: ${process.stderr}" },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Configuration Issues

### Profile Problems

**Problem**: Wrong profile used or profile not found

```
ERROR: Unknown profile 'test'
```

**Solutions**:

1. **Check available profiles**:

   ```bash
   nextflow config -show-profiles
   ```

2. **Verify profile definition**:

   ```groovy
   // In nextflow.config
   profiles {
       test {
           params {
               max_memory = '6.GB'
               max_cpus = 2
           }
       }
   }
   ```

3. **Use explicit configuration**:
   ```bash
   nf-test test -c tests/test.config
   ```

### Parameter Conflicts

**Problem**: Parameter values not being applied

**Solutions**:

1. **Check parameter precedence**:

   ```bash
   # Parameters are applied in this order:
   # 1. Command line: --param value
   # 2. Config files: params.param = value
   # 3. Default values in script
   ```

2. **Debug parameter values**:
   ```groovy
   test("Parameter debug") {
       when {
           params {
               custom_param = 'test_value'
           }
           process {
               """
               echo "Parameter value: ${params.custom_param}"
               echo "All params: ${params}"

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
               { assert process.stdout.contains('test_value') }
           )
       }
   }
   ```

## Environment Issues

### Java/Groovy Problems

**Problem**: Groovy compilation errors

```
ERROR: Script compilation error
groovy.lang.MissingMethodException: No signature of method...
```

**Solutions**:

1. **Check Groovy syntax**:

   ```groovy
   // Correct closure syntax
   { assert process.success }

   // Not: assert process.success (missing closure braces)
   ```

2. **Validate method calls**:

   ```groovy
   // Check available methods
   println process.class.methods*.name.sort()
   ```

3. **Java version compatibility**:
   ```bash
   java -version
   echo $JAVA_HOME
   ```

### Path and Permissions

**Problem**: Permission denied or path issues

```
ERROR: Cannot read file: permission denied
```

**Solutions**:

1. **Check file permissions**:

   ```bash
   ls -la tests/data/
   chmod 644 tests/data/*.fastq.gz
   ```

2. **Verify working directory**:
   ```groovy
   test("Path debug") {
       when {
           process {
               """
               echo "Working directory: \$(pwd)"
               echo "File listing:"
               ls -la

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

## Performance Issues

### Slow Tests

**Problem**: Tests take too long to complete

**Solutions**:

1. **Use smaller test data**:

   ```groovy
   // Create minimal datasets
   setup {
       def quickFile = new File("tests/quick_test.fastq")
       quickFile.text = (1..100).collect { i ->
           "@read${i}\nATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIII"
       }.join('\n')
   }
   ```

2. **Optimize resource allocation**:

   ```groovy
   // In test config
   process {
       withName: 'SLOW_PROCESS' {
           cpus = 4
           memory = '8.GB'
       }
   }
   ```

3. **Use caching**:
   ```bash
   # Enable Nextflow caching
   nf-test test -resume
   ```

### Memory Leaks

**Problem**: Memory usage increases over time

**Solutions**:

1. **Monitor memory usage**:

   ```bash
   # Run tests with memory monitoring
   while true; do
       echo "$(date): $(free -h | grep Mem)"
       sleep 10
   done &

   nf-test test
   ```

2. **Clean up between tests**:
   ```groovy
   cleanup {
       // Clean temporary files
       new File('temp_files').deleteDir()

       // Force garbage collection
       System.gc()
   }
   ```

## Data Issues

### Corrupted Test Data

**Problem**: Test data appears corrupted

```
ERROR: Invalid FASTQ format
```

**Solutions**:

1. **Validate test data integrity**:

   ```bash
   # Check file integrity
   md5sum tests/data/*.fastq.gz

   # Validate FASTQ format
   zcat test.fastq.gz | head -4
   ```

2. **Re-download test data**:

   ```bash
   # Clean and re-download
   rm -rf tests/data/
   curl -L https://github.com/nf-core/test-datasets/archive/modules.tar.gz | \
       tar -xz --strip-components=1 -C tests/data/
   ```

3. **Generate fresh test data**:
   ```groovy
   setup {
       // Generate valid test data
       def validFastq = new File("tests/generated.fastq")
       validFastq.text = """@read1
   ATCGATCGATCGATCGATCGATCGATCGATCG
   +
   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
   @read2
   GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
   +
   JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
   """
   }
   ```

### Missing Dependencies

**Problem**: Required tools not available

```
ERROR: Command not found: samtools
```

**Solutions**:

1. **Check tool availability**:

   ```groovy
   test("Tool availability check") {
       when {
           process {
               """
               # Check for required tools
               for tool in samtools bcftools bwa; do
                   if ! command -v \$tool >/dev/null 2>&1; then
                       echo "ERROR: \$tool not found in PATH"
                       exit 1
                   fi
                   echo "\$tool: \$(\$tool --version | head -1)"
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
               { assert process.stdout.contains('samtools:') }
           )
       }
   }
   ```

2. **Use container fallback**:
   ```groovy
   // In nextflow.config
   process {
       withName: 'TOOL' {
           container = 'biocontainers/samtools:1.15.1--h1170115_0'
       }
   }
   ```

## Assertion Failures

### Complex Assertion Debugging

**Problem**: Assertion fails but reason unclear

```
Assertion failed: assert process.out.bam.size() == 2
```

**Solutions**:

1. **Add detailed assertions**:

   ```groovy
   then {
       assertAll(
           { assert process.success },
           {
               println "Number of BAM outputs: ${process.out.bam.size()}"
               println "BAM outputs: ${process.out.bam}"
               process.out.bam.eachWithIndex { item, index ->
                   println "Output ${index}: ${item}"
               }
               assert process.out.bam.size() == 2
           }
       )
   }
   ```

2. **Use descriptive assertion messages**:

   ```groovy
   {
       def actualSize = process.out.bam.size()
       def expectedSize = 2
       assert actualSize == expectedSize,
              "Expected ${expectedSize} BAM files, but got ${actualSize}. Outputs: ${process.out.bam}"
   }
   ```

3. **Check intermediate outputs**:
   ```groovy
   {
       // Verify each output channel
       assert process.out.containsKey('bam'), "Missing 'bam' output channel"
       assert process.out.bam instanceof List, "BAM output is not a list"
       assert !process.out.bam.isEmpty(), "BAM output list is empty"

       process.out.bam.each { meta, bam ->
           assert meta instanceof Map, "Meta is not a Map: ${meta}"
           assert meta.containsKey('id'), "Meta missing 'id' field: ${meta}"
           assert bam.exists(), "BAM file does not exist: ${bam}"
       }
   }
   ```

## CI/CD Issues

### GitHub Actions Failures

**Problem**: Tests pass locally but fail in CI

**Solutions**:

1. **Check environment differences**:

   ```yaml
   - name: Debug environment
     run: |
       echo "OS: $(uname -a)"
       echo "Java: $(java -version)"
       echo "Available memory: $(free -h)"
       echo "Available disk: $(df -h)"
       echo "CPU info: $(nproc)"
   ```

2. **Use consistent test data**:

   ```yaml
   - name: Cache test data
     uses: actions/cache@v3
     with:
       path: tests/data/
       key: test-data-${{ hashFiles('tests/data.md5') }}
   ```

3. **Add CI-specific configuration**:
   ```groovy
   // In nextflow.config
   profiles {
       github_actions {
           params {
               max_memory = '6.GB'
               max_cpus = 2
               max_time = '1.h'
           }
           process.container = 'ubuntu:20.04'
       }
   }
   ```

### Docker Issues in CI

**Problem**: Docker commands fail in CI

**Solutions**:

1. **Check Docker daemon**:

   ```yaml
   - name: Check Docker
     run: |
       docker --version
       docker info
       docker ps
   ```

2. **Use pre-built images**:
   ```yaml
   - name: Pull required images
     run: |
       docker pull biocontainers/fastqc:0.11.9--0
       docker pull biocontainers/samtools:1.15.1--h1170115_0
   ```

## General Debugging Tips

### Enable Debug Mode

```bash
# Enable all debugging options
export NXF_DEBUG=1
export NEXTFLOW_DEBUG=1

# Run with maximum verbosity
nf-test test --verbose --debug -with-trace -with-report
```

### Use Test Isolation

```groovy
// Isolate problematic tests
test("Isolated test") {
    cleanup {
        // Clean up after test
        new File('work').deleteDir()
        new File('.nextflow').deleteDir()
    }

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

### Create Minimal Reproduction Cases

```groovy
// Simplify test to isolate issue
test("Minimal reproduction") {
    when {
        process {
            """
            # Simplest possible test case
            echo "test" > output.txt

            input[0] = [
                [ id:'test' ],
                []
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert new File('output.txt').exists() }
        )
    }
}
```

## Getting Help

### Community Resources

1. **nf-core Slack**: `#help` channel
2. **GitHub Issues**: nf-test repository issues
3. **nf-core Documentation**: https://nf-co.re/docs
4. **Nextflow Documentation**: https://nextflow.io/docs

### Reporting Issues

When reporting issues, include:

1. **nf-test version**: `nf-test --version`
2. **Nextflow version**: `nextflow -version`
3. **System information**: OS, Java version
4. **Complete error message**
5. **Minimal reproduction case**
6. **Configuration files used**

### Example Issue Report

```markdown
## Bug Report

**nf-test version**: 0.8.4
**Nextflow version**: 23.10.0
**OS**: Ubuntu 20.04
**Java**: OpenJDK 17.0.7

### Error Message
```

ERROR: Test failed with snapshot mismatch
Expected: test.html:md5:abc123
Actual: test.html:md5:def456

```

### Steps to Reproduce
1. Run `nf-test test modules/fastqc/tests/main.nf.test`
2. Error occurs on "Single-end reads" test

### Configuration
- Profile: docker
- Test data: nf-core modules test datasets

### Expected Behavior
Test should pass with matching snapshots

### Additional Context
- Test passes locally but fails in CI
- Only affects FastQC module
```

This comprehensive troubleshooting guide should help you debug most common issues encountered when testing with nf-test in nf-core environments.
