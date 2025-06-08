---
title: "12. External Commands and Tools"
subtitle: Testing with external dependencies
weight: 120
---

## Testing with External Tools

### Tool Availability Checks

```groovy
test("Tool availability") {
    when {
        process {
            """
            // Check if required tools are available
            which samtools || exit 1
            which bcftools || exit 1
            
            input[0] = [
                [ id:'test' ],
                file('test.bam', checkIfExists: true)
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

### Version-Specific Testing

```groovy
test("Tool version compatibility") {
    when {
        process {
            """
            // Check tool version
            SAMTOOLS_VERSION=\$(samtools --version | head -n1 | cut -d' ' -f2)
            echo "Using samtools version: \$SAMTOOLS_VERSION"
            
            # Version-specific logic
            if [[ "\$SAMTOOLS_VERSION" > "1.10" ]]; then
                SORT_ARGS="--threads ${task.cpus}"
            else
                SORT_ARGS="-@ ${task.cpus}"
            fi
            
            input[0] = [
                [ id:'test', samtools_version: "\$SAMTOOLS_VERSION" ],
                file('test.bam', checkIfExists: true)
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

## Container Testing

### Multi-Container Testing

```groovy
nextflow_process {
    name "Test Process with Multiple Containers"
    script "../main.nf"
    process "MULTI_TOOL"
    
    test("Docker container") {
        when {
            process {
                """
                container 'biocontainers/samtools:1.15.1--h1170115_0'
                
                input[0] = [
                    [ id:'test' ],
                    file('test.bam', checkIfExists: true)
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
    
    test("Singularity container") {
        when {
            process {
                """
                container 'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'
                
                input[0] = [
                    [ id:'test' ],
                    file('test.bam', checkIfExists: true)
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
}
```

### Container Environment Testing

```groovy
test("Container environment validation") {
    when {
        process {
            """
            // Validate container environment
            echo "Container info:"
            uname -a
            echo "Available tools:"
            which samtools && samtools --version
            which bcftools && bcftools --version
            echo "Environment variables:"
            env | grep -E "(PATH|LD_LIBRARY_PATH)" | sort
            
            input[0] = [
                [ id:'test' ],
                file('test.bam', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.stdout.contains("samtools") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Conda Environment Testing

### Environment Validation

```groovy
test("Conda environment") {
    when {
        process {
            """
            // Activate conda environment
            source activate nf-core-module-env
            
            // Validate environment
            conda list | grep samtools
            python --version
            
            input[0] = [
                [ id:'test' ],
                file('test.bam', checkIfExists: true)
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

### Package Version Testing

```groovy
test("Package versions") {
    when {
        process {
            """
            // Check specific package versions
            python -c "import pandas; print(f'pandas: {pandas.__version__}')"
            python -c "import numpy; print(f'numpy: {numpy.__version__}')"
            R --version | head -n1
            
            input[0] = [
                [ id:'test' ],
                file('data.csv', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.stdout.contains("pandas:") },
            { assert process.stdout.contains("numpy:") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Database and Reference Testing

### Database Availability

```groovy
test("Database access") {
    when {
        process {
            """
            // Check database availability
            if [[ -n "\${NCBI_API_KEY:-}" ]]; then
                echo "NCBI API key available"
                # Test API access
                curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi" | head
            else
                echo "No NCBI API key, using local database"
            fi
            
            # Check local database
            if [[ -f "\${params.blast_db}/nt.nal" ]]; then
                echo "Local BLAST database found"
            else
                echo "Local BLAST database not found, downloading..."
                # Download minimal test database
            fi
            
            input[0] = [
                [ id:'test' ],
                file('query.fasta', checkIfExists: true)
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

### Reference Genome Testing

```groovy
test("Reference genome validation") {
    when {
        process {
            """
            // Validate reference genome
            REFERENCE="${params.genome_fasta}"
            
            if [[ ! -f "\$REFERENCE" ]]; then
                echo "ERROR: Reference genome not found: \$REFERENCE"
                exit 1
            fi
            
            # Check reference format
            if ! head -n1 "\$REFERENCE" | grep -q "^>"; then
                echo "ERROR: Invalid FASTA format"
                exit 1
            fi
            
            # Check if indexed
            if [[ ! -f "\${REFERENCE}.fai" ]]; then
                echo "Creating FASTA index"
                samtools faidx "\$REFERENCE"
            fi
            
            input[0] = [
                [ id:'test' ],
                file('reads.fastq.gz', checkIfExists: true)
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

## Network and API Testing

### API Connectivity

```groovy
test("External API access") {
    when {
        process {
            """
            // Test API connectivity
            API_URL="https://rest.ensembl.org"
            
            # Check if API is accessible
            if curl -s --connect-timeout 10 "\$API_URL/info/ping" | grep -q "ping"; then
                echo "Ensembl API accessible"
                API_AVAILABLE=true
            else
                echo "Ensembl API not accessible, using cached data"
                API_AVAILABLE=false
            fi
            
            input[0] = [
                [ id:'test', api_available: \$API_AVAILABLE ],
                file('gene_ids.txt', checkIfExists: true)
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

### Download Testing

```groovy
test("External data download") {
    when {
        process {
            """
            // Test data download with fallback
            DOWNLOAD_URL="${params.reference_url}"
            LOCAL_CACHE="${params.cache_dir}/reference.fa"
            
            if [[ -f "\$LOCAL_CACHE" ]]; then
                echo "Using cached reference"
                cp "\$LOCAL_CACHE" reference.fa
            else
                echo "Downloading reference from \$DOWNLOAD_URL"
                if curl -s -L -o reference.fa "\$DOWNLOAD_URL"; then
                    echo "Download successful"
                    # Cache for future use
                    mkdir -p "\${params.cache_dir}"
                    cp reference.fa "\$LOCAL_CACHE"
                else
                    echo "Download failed, using minimal test reference"
                    echo ">test_seq" > reference.fa
                    echo "ATCGATCGATCGATCG" >> reference.fa
                fi
            fi
            
            input[0] = [
                [ id:'test' ],
                file('reads.fastq.gz', checkIfExists: true)
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

## License and Compliance Testing

### License Validation

```groovy
test("Software licensing") {
    when {
        process {
            """
            // Check software licenses
            echo "Checking software licenses..."
            
            # Check if commercial software is properly licensed
            if command -v matlab >/dev/null 2>&1; then
                echo "MATLAB found, checking license..."
                matlab -nodisplay -nodesktop -r "license('checkout','statistics_toolbox'); exit" 2>&1 | tee matlab_license.log
                if grep -q "checked out" matlab_license.log; then
                    echo "MATLAB license valid"
                    MATLAB_AVAILABLE=true
                else
                    echo "MATLAB license invalid, using alternative"
                    MATLAB_AVAILABLE=false
                fi
            else
                echo "MATLAB not found, using open-source alternative"
                MATLAB_AVAILABLE=false
            fi
            
            input[0] = [
                [ id:'test', matlab_available: \$MATLAB_AVAILABLE ],
                file('data.mat', checkIfExists: true)
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

## Performance and Resource Testing

### Resource Monitoring

```groovy
test("External tool resource usage") {
    when {
        process {
            """
            # Monitor resource usage of external tools
            echo "Starting resource monitoring..."
            
            # Start background monitoring
            (while true; do
                ps aux | grep -E "(samtools|bcftools)" | grep -v grep
                free -h
                sleep 5
            done) > resource_monitor.log &
            MONITOR_PID=\$!
            
            # Run the actual process
            samtools sort -@ ${task.cpus} input.bam -o output.bam
            
            # Stop monitoring
            kill \$MONITOR_PID 2>/dev/null || true
            
            # Analyze resource usage
            echo "Peak memory usage:"
            grep -E "samtools.*sort" resource_monitor.log | awk '{print \$6}' | sort -n | tail -1
            
            input[0] = [
                [ id:'test' ],
                file('input.bam', checkIfExists: true)
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

### Timeout and Retry Logic

```groovy
test("External tool with timeout") {
    when {
        process {
            """
            # Function to run command with timeout
            run_with_timeout() {
                local cmd="\$1"
                local timeout="\$2"
                
                timeout "\$timeout" bash -c "\$cmd" || {
                    echo "Command timed out after \$timeout seconds"
                    return 1
                }
            }
            
            # Function to retry command
            retry_command() {
                local cmd="\$1"
                local max_attempts=3
                local attempt=1
                
                while [[ \$attempt -le \$max_attempts ]]; do
                    echo "Attempt \$attempt of \$max_attempts"
                    if run_with_timeout "\$cmd" "300s"; then
                        echo "Command succeeded on attempt \$attempt"
                        return 0
                    else
                        echo "Command failed on attempt \$attempt"
                        attempt=\$((attempt + 1))
                        sleep 10
                    fi
                done
                
                echo "Command failed after \$max_attempts attempts"
                return 1
            }
            
            # Use retry logic for unreliable external command
            retry_command "external_unreliable_tool --input input.txt --output output.txt"
            
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

## Mock and Stub Testing

### Mocking External Services

```groovy
test("Mock external service") {
    setup {
        // Create mock service response
        def mockResponse = new File("tests/mock_api_response.json")
        mockResponse.text = '''
        {
            "status": "success",
            "data": {
                "gene_id": "ENSG00000123456",
                "gene_name": "TEST_GENE",
                "chromosome": "1"
            }
        }
        '''
    }
    
    when {
        process {
            """
            # Mock external API call
            if [[ "\${params.use_mock_api}" == "true" ]]; then
                echo "Using mock API response"
                cp tests/mock_api_response.json api_response.json
            else
                echo "Calling real API"
                curl -s "https://api.example.com/gene/ENSG00000123456" > api_response.json
            fi
            
            input[0] = [
                [ id:'test' ],
                file('gene_list.txt', checkIfExists: true)
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

### Tool Stubbing

```groovy
test("Stubbed external tool") {
    setup {
        // Create stub for external tool
        def stubScript = new File("tests/stubs/external_tool")
        stubScript.text = '''#!/bin/bash
echo "STUB: external_tool called with arguments: $@"
echo "Creating mock output..."
echo "mock_result" > $2  # Assuming $2 is output file
exit 0
'''
        stubScript.setExecutable(true)
    }
    
    when {
        process {
            """
            # Add stub to PATH
            export PATH="tests/stubs:\$PATH"
            
            # Verify stub is used
            which external_tool
            
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
            { assert process.stdout.contains("STUB:") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Error Handling and Fallbacks

### Graceful Degradation

```groovy
test("Graceful tool degradation") {
    when {
        process {
            """
            # Try to use preferred tool, fallback to alternative
            if command -v fast_tool >/dev/null 2>&1; then
                echo "Using fast_tool"
                fast_tool --input input.txt --output output.txt
            elif command -v slow_tool >/dev/null 2>&1; then
                echo "fast_tool not available, using slow_tool"
                slow_tool --input input.txt --output output.txt
            else
                echo "No suitable tool found, using built-in method"
                # Built-in fallback implementation
                cat input.txt | sed 's/pattern/replacement/g' > output.txt
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

## Next Steps

Continue to [Custom Utility Classes](./13_custom_utility_classes.md) to learn about creating reusable test utilities. 