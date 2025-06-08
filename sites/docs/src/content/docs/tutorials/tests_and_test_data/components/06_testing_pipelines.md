---
title: "6. Testing Pipelines"
subtitle: End-to-end pipeline testing
weight: 60
---

## Pipeline Test Structure

Pipeline tests verify end-to-end functionality with complete datasets. Here's the standard structure used in nf-core pipelines:

```groovy
nextflow_pipeline {
    name "Test pipeline"
    script "../main.nf"
    config "./nextflow.config"
    tag "pipeline"
    tag "cpu"

    test("-profile test") {
        when {
            params {
                outdir = "$outputDir"
            }
        }
        then {
            // stable_name: All files + folders in ${params.outdir}/ with a stable name
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
            // stable_path: All files in ${params.outdir}/ with stable content
            def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
            // binary_files: All binary files that need checksum validation
            def binary_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam', '**/*.bai', '**/*.vcf.gz'])
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    // Number of tasks
                    workflow.trace.succeeded().size(),
                    // pipeline versions.yml file for multiqc from which Nextflow version is removed because we test pipelines on multiple Nextflow versions
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                    // All stable path names
                    stable_name,
                    // All files with stable contents
                    stable_path,
                    // All binary files with checksums
                    binary_files.collect{ file -> [ file.getName(), file.getName().endsWith('.bam') ? bam(file.toString()).getReadsMD5() : file.md5() ] }
                ).match() }
            )
        }
    }
}
```

## Mandatory Pipeline Tests

### Required default.nf.test

Every nf-core pipeline **MUST** have a `default.nf.test` file in the pipeline root. Here's the standard pattern:

```groovy
// default.nf.test
nextflow_pipeline {
    name "Test pipeline"
    script "../main.nf"
    config "./nextflow.config"
    tag "pipeline"
    tag "cpu"

    test("-profile test") {
        when {
            params {
                outdir = "$outputDir"
            }
        }
        then {
            // Standard nf-core pattern with comprehensive snapshot coverage
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
            def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
            def binary_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam', '**/*.bai', '**/*.vcf.gz'])
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.trace.succeeded().size(),
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                    stable_name,
                    stable_path,
                    binary_files.collect{ file -> [ file.getName(), file.getName().endsWith('.bam') ? bam(file.toString()).getReadsMD5() : file.md5() ] }
                ).match() }
            )
        }
    }
}
```

## Basic Pipeline Testing

### Minimal Test Configuration

```groovy
nextflow_pipeline {
    name "Test Pipeline"
    script "main.nf"
    
    test("Basic pipeline run") {
        when {
            params {
                input = 'tests/minimal_samplesheet.csv'
                outdir = './test_results'
                genome = 'test'
            }
        }
        then {
            assertAll(
                { assert workflow.success },
                { assert new File("$outputDir/pipeline_info/pipeline_dag.html").exists() },
                { assert new File("$outputDir/multiqc/multiqc_report.html").exists() }
            )
        }
    }
}
```

### Using Test Profiles

```groovy
nextflow_pipeline {
    name "Test Pipeline with Profile"
    script "main.nf"
    profile "test"
    
    test("Test profile execution") {
        when {
            // Profile provides default test parameters
        }
        then {
            assertAll(
                { assert workflow.success },
                { assert workflow.trace.succeeded().size() > 0 }
            )
        }
    }
}
```

## Metro Map Path Testing

### Testing Each Pipeline Path

Test each major pathway through your pipeline's "metro map". Here are examples for common pipeline branches:

```groovy
// Test main pipeline path (default)
test("-profile test") {
    when {
        params {
            outdir = "$outputDir"
        }
    }
    then {
        def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
        def binary_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam', '**/*.bai', '**/*.vcf.gz'])
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.trace.succeeded().size(),
                removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                stable_name,
                stable_path,
                binary_files.collect{ file -> [ file.getName(), file.getName().endsWith('.bam') ? bam(file.toString()).getReadsMD5() : file.md5() ] }
            ).match() }
        )
    }
}

// Test alternate aligner path - bwameth
test("Params: bwameth") {
    when {
        params {
            aligner = "bwameth"
            outdir = "$outputDir"
        }
    }
    then {
        def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
        def bam_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam'])
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.trace.succeeded().size(),
                removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                stable_name,
                stable_path,
                bam_files.collect{ file -> [ file.getName(), bam(file.toString()).getReadsMD5() ] }
            ).match() }
        )
    }
}

// Test targeted sequencing path
test("Params: bismark | run_targeted_sequencing") {
    when {
        params {
            outdir = "$outputDir"
            run_targeted_sequencing = true
            target_regions_file = "https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome_target_regions.bed"
        }
    }
    then {
        def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
        def bam_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam'])
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.trace.succeeded().size(),
                removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                stable_name,
                stable_path,
                bam_files.collect{ file -> [ file.getName(), bam(file.toString()).getReadsMD5() ] }
            ).match() }
        )
    }
}

// Test skip options path
test("Params: bismark | skip_trimming") {
    when {
        params {
            skip_trimming = true
            outdir = "$outputDir"
        }
    }
    then {
        def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
        def bam_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam'])
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.trace.succeeded().size(),
                removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                stable_name,
                stable_path,
                bam_files.collect{ file -> [ file.getName(), bam(file.toString()).getReadsMD5() ] }
            ).match() }
        )
    }
}
```

## Testing Different Pipeline Configurations

### Testing Multiple Profiles

```groovy
test("Test profile - minimal") {
    profile "test"
    
    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(workflow.out).match("test_profile") }
        )
    }
}

test("Full test profile") {
    profile "test_full"
    
    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(workflow.out).match("test_full_profile") }
        )
    }
}
```

### Crucial Parameter Testing

Test the most important parameters that significantly affect pipeline behavior:

```groovy
// Test critical boolean parameters
test("Skip QC steps") {
    when {
        params {
            input = 'tests/samplesheet.csv'
            skip_fastqc = true
            skip_multiqc = true
            outdir = './results_no_qc'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert !new File("$outputDir/multiqc/").exists() },
            { assert workflow.trace.succeeded().size() > 0 }
        )
    }
}

// Test alternative tool selection
test("Alternative aligner") {
    when {
        params {
            input = 'tests/samplesheet.csv'
            aligner = 'star'
            outdir = './results_star'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(workflow.out).match("star_aligner") },
            { assert workflow.trace.succeeded().size() > 0 }
        )
    }
}

// Test reference genome parameter
test("Different genome build") {
    when {
        params {
            input = 'tests/samplesheet.csv'
            genome = 'GRCh37'
            outdir = './results_grch37'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(workflow.out).match("grch37_genome") },
            { assert workflow.trace.succeeded().size() > 0 }
        )
    }
}
```

## Using nft-utils Plugin

### Installing and Configuring nft-utils

Add to your `nf-test.config`:

```groovy
config {
    plugins {
        load "nft-utils@0.0.4"
    }
    
    testsDir "tests"
    workDir ".nf-test"
    configFile "tests/nextflow.config"
}
```

### Essential Pipeline-Level Snapshot Pattern

This is the **recommended pattern** for comprehensive pipeline testing:

```groovy
nextflow_pipeline {
    name "Test Pipeline Default"
    script "main.nf"
    tag "pipeline"
    tag "default"

    test("Pipeline with stable snapshots") {
        when {
            params {
                input = 'tests/samplesheet.csv'
                outdir = "$outputDir"
                genome = 'GRCh37'
            }
        }
        then {
            // stable_name: All files + folders in ${params.outdir}/ with stable names
            def stable_name = getAllFilesFromDir(
                params.outdir, 
                relative: true, 
                includeDir: true, 
                ignore: ['pipeline_info/*.{html,json,txt}', 'pipeline_info/execution_*.{html,txt}']
            )
            
            // stable_path: All files in ${params.outdir}/ with stable content
            def stable_path = getAllFilesFromDir(
                params.outdir, 
                ignoreFile: 'tests/.nftignore'
            )
            
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    // Number of successful tasks
                    workflow.trace.succeeded().size(),
                    // Pipeline versions.yml file with Nextflow version removed
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                    // All stable file/folder names (relative paths)
                    stable_name,
                    // All files with stable contents
                    stable_path
                ).match() }
            )
        }
    }
}
```

### Setting up .nftignore for Unstable Files

Create `tests/.nftignore` to exclude files with unstable content but stable names. Here's an example:

```gitignore
.DS_Store
*/alignments/logs/*.txt
*/{alignments,deduplicated}/*.{bam,bam.bai}
*/deduplicated/logs/*.txt
*/{reports,summary}/*.{html,txt}
*/deduplicated/picard_metrics/*.txt
fastqc/*.html
fastqc/zips/*.zip
multiqc/*/multiqc_data/*.{log,json}
multiqc/*/multiqc_data/multiqc_fastqc.txt
multiqc/*/multiqc_data/multiqc_general_stats.txt
multiqc/*/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt
multiqc/*/multiqc_data/multiqc_software_versions.txt
multiqc/*/multiqc_data/multiqc_sources.txt
multiqc/*/multiqc_report.html
multiqc/*/multiqc_plots/{pdf,png,svg}/*.{pdf,png,svg}
pipeline_info/*.{html,json,txt,yml}
*/qualimap/bamqc/*/qualimapReport.html
*/qualimap/bamqc/**/*.{pdf,png,svg}
*/qualimap/bamqc/*/css/*
qualimap/**/images_qualimapReport/*
qualimap/**/raw_data_qualimapReport/*
trimgalore/fastqc/*.html
trimgalore/fastqc/zips/*.zip
trimgalore/logs/*.txt
enrichment_metrics/*
```

### Advanced nft-utils Usage

#### Removing Specific YAML Keys

Remove specific entries from YAML files that vary between runs:

```groovy
test("Pipeline without version dependencies") {
    when {
        params {
            input = 'tests/samplesheet.csv'
            outdir = "$outputDir"
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                // Remove entire Workflow section
                removeFromYamlMap("$outputDir/pipeline_info/*_versions.yml", "Workflow"),
                // Or remove specific subkey
                removeFromYamlMap("$outputDir/pipeline_info/*_versions.yml", "Workflow", "Nextflow"),
                workflow.trace.succeeded().size()
            ).match() }
        )
    }
}
```

#### Named Parameters for getAllFilesFromDir

The function supports flexible named parameters:

```groovy
// Example with all named parameters
def comprehensive_capture = getAllFilesFromDir(
    params.outdir,
    includeDir: true,          // Include directory names
    relative: true,            // Use relative paths
    ignore: [                  // Patterns to exclude
        'pipeline_info/*.{html,json,txt}',
        '**/*.log',
        '**/work/**'
    ],
    ignoreFile: 'tests/.nftignore',  // Additional ignore patterns from file
    include: ['*', '**/*']     // Patterns to include (default: all)
)
```

#### Wildcard Pattern Support

Both `removeNextflowVersion()` and `removeFromYamlMap()` support wildcards:

```groovy
// Process all version files in any subdirectory
removeNextflowVersion("$outputDir/**/versions.yml")

// Remove from multiple files
removeFromYamlMap("$outputDir/pipeline_info/*_mqc_versions.yml", "Workflow", "Nextflow")
```

### Real-World Pipeline Test Examples

#### Complete Example Test

```groovy
nextflow_pipeline {
    name "Test Example Pipeline"
    script "main.nf"
    profile "test"
    tag "pipeline"

    test("Test profile execution") {
        when {
            params {
                outdir = "$outputDir"
            }
        }
        then {
            // Separate stable names and stable content
            def stable_name = getAllFilesFromDir(
                params.outdir, 
                relative: true, 
                includeDir: true, 
                ignore: [
                    'pipeline_info/*.{html,json,txt}',
                    'pipeline_info/execution_*.{html,txt}',
                    '**/fastqc/**/*.html',
                    '**/multiqc/**/*.html'
                ]
            )
            
            def stable_path = getAllFilesFromDir(
                params.outdir, 
                ignoreFile: 'tests/.nftignore'
            )
            
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.trace.succeeded().size(),
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                    stable_name,
                    stable_path
                ).match() }
            )
        }
    }
}
```

#### Multi-path Pipeline Testing

```groovy
nextflow_pipeline {
    name "Test Multiple Pipeline Paths"
    script "main.nf"
    
    test("Bismark aligner path") {
        when {
            params {
                input = 'tests/samplesheet.csv'
                outdir = "$outputDir"
                aligner = 'bismark'
            }
        }
        then {
            def stable_outputs = getAllFilesFromDir(
                params.outdir,
                relative: true,
                ignore: ['pipeline_info/execution_*.{html,txt}']
            )
            
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.trace.succeeded().size(),
                    removeNextflowVersion("$outputDir/pipeline_info/*_versions.yml"),
                    stable_outputs
                ).match("bismark_path") }
            )
        }
    }
    
    test("BWA-meth aligner path") {
        when {
            params {
                input = 'tests/samplesheet.csv'
                outdir = "$outputDir"
                aligner = 'bwameth'
            }
        }
        then {
            def stable_outputs = getAllFilesFromDir(
                params.outdir,
                relative: true,
                ignore: ['pipeline_info/execution_*.{html,txt}']
            )
            
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.trace.succeeded().size(),
                    removeNextflowVersion("$outputDir/pipeline_info/*_versions.yml"),
                    stable_outputs
                ).match("bwameth_path") }
            )
        }
    }
}
```

### Best Practices for nft-utils

1. **Always set outdir to `$outputDir`** - ensures consistent test directory structure
2. **Use meaningful snapshot tags** - helps identify specific test scenarios
3. **Separate stable names from stable content** - provides better debugging when tests fail
4. **Include task count verification** - `workflow.trace.succeeded().size()` ensures execution completeness
5. **Remove version dependencies** - use `removeNextflowVersion()` for reproducible tests
6. **Maintain .nftignore files** - keep them updated as pipeline outputs evolve

### Troubleshooting nft-utils

#### When snapshots fail to match:

1. **Check .nftignore patterns** - ensure unstable files are properly excluded
2. **Verify wildcard patterns** - make sure they match your actual file structure
3. **Review relative vs absolute paths** - use `relative: true` for portable tests
4. **Update ignore patterns** - add new unstable files discovered during testing

#### Example debugging approach:

```groovy
// Temporary debugging - capture everything to see what's unstable
def debug_all = getAllFilesFromDir(params.outdir, relative: true, includeDir: true)
println "All files found: ${debug_all}"

// Then refine ignore patterns based on output
```

## Advanced Pipeline Testing

### Testing with Custom Input

```groovy
test("Custom samplesheet validation") {
    setup {
        // Create custom samplesheet
        def samplesheet = new File("tests/custom_samples.csv")
        samplesheet.text = """sample,fastq_1,fastq_2,strandedness
sample1,test_1.fastq.gz,test_2.fastq.gz,auto
sample2,test2_1.fastq.gz,test2_2.fastq.gz,forward
"""
    }
    
    when {
        params {
            input = 'tests/custom_samples.csv'
            outdir = './custom_results'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert workflow.trace.succeeded().size() >= 4 } // At least 4 processes per sample
        )
    }
}
```

### Resource Testing

```groovy
test("Resource limits") {
    when {
        params {
            input = 'tests/samplesheet.csv'
            max_cpus = 2
            max_memory = '4.GB'
            max_time = '2.h'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            // Verify no process exceeded resource limits
            { 
                workflow.trace.each { task ->
                    assert task.get('cpus') <= 2
                    assert task.get('memory') <= 4_000_000_000 // 4GB in bytes
                }
            }
        )
    }
}
```

## Testing Pipeline Outputs

### Comprehensive Output Verification

```groovy
then {
    assertAll(
        { assert workflow.success },
        
        // Check main outputs exist
        { assert new File("$outputDir/star/").exists() },
        { assert new File("$outputDir/salmon/").exists() },
        { assert new File("$outputDir/multiqc/").exists() },
        
        // Check specific files
        { assert new File("$outputDir/multiqc/multiqc_report.html").exists() },
        { assert new File("$outputDir/pipeline_info/execution_report.html").exists() },
        
        // Verify output structure
        { 
            def samples = ['sample1', 'sample2']
            samples.each { sample ->
                assert new File("$outputDir/star/${sample}/").exists()
                assert new File("$outputDir/salmon/${sample}/").exists()
            }
        },
        
        // Snapshot key outputs
        { assert snapshot(
            workflow.out.multiqc_report,
            workflow.out.software_versions
        ).match() }
    )
}
```

### Selective Output Testing

```groovy
then {
    assertAll(
        { assert workflow.success },
        
        // Test specific output channels
        { assert workflow.out.star_aligned.size() == 2 },
        { assert workflow.out.salmon_results.size() == 2 },
        
        // Verify file extensions
        {
            workflow.out.star_aligned.each { meta, bam ->
                assert bam.toString().endsWith('.bam')
                assert meta.id in ['sample1', 'sample2']
            }
        }
    )
}
```

## Stub Testing for Pipelines

```groovy
test("Pipeline stub test") {
    options "-stub-run"
    
    when {
        params {
            input = 'tests/samplesheet_large.csv'
            outdir = './stub_results'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            { assert workflow.trace.succeeded().size() > 0 },
            // Verify stub outputs were created
            { assert snapshot(workflow.out).match("stub_outputs") }
        )
    }
}
```

## Error Testing

### Invalid Input Testing

```groovy
test("Invalid samplesheet") {
    when {
        params {
            input = 'tests/invalid_samplesheet.csv'
            outdir = './error_results'
        }
    }
    then {
        assertAll(
            { assert workflow.failed },
            { assert workflow.stderr.contains("ERROR: Validation of pipeline parameters failed!") }
        )
    }
}
```

### Missing Parameter Testing

```groovy
test("Missing required parameter") {
    when {
        params {
            outdir = './results'
            // Missing required 'input' parameter
        }
    }
    then {
        assertAll(
            { assert workflow.failed },
            { assert workflow.stderr.contains("Missing required parameter") }
        )
    }
}
```

## Performance Testing

### Execution Time Verification

```groovy
test("Performance benchmark") {
    when {
        params {
            input = 'tests/benchmark_samplesheet.csv'
            outdir = './benchmark_results'
        }
    }
    then {
        assertAll(
            { assert workflow.success },
            // Verify execution completed within reasonable time
            { assert workflow.duration.toMillis() < 600_000 }, // 10 minutes
            { assert workflow.trace.succeeded().size() >= 10 }
        )
    }
}
```

## Next Steps

Continue to [nf-test Assertions](./07_assertions.md) to learn about comprehensive assertion patterns and verification techniques. 

### Setting up nf-test Configuration

Add to your `nf-test.config` based on nf-core standard implementation:

```groovy
config {
    // location for all nf-test tests
    testsDir "."

    // nf-test directory including temporary files for each test
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"

    // location of an optional nextflow.config file specific for executing tests
    configFile "tests/nextflow.config"

    // ignore tests coming from the nf-core/modules repo
    ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'

    // run all test with defined profile(s) from the main nextflow.config
    profile "test"

    // list of filenames or patterns that should be trigger a full test run
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'tests/nextflow.config', 'tests/.nftignore'

    // load the necessary plugins
    plugins {
        load "nft-bam@0.5.0"
        load "nft-utils@0.0.3"
    }
}
```

### Test Configuration

Create `tests/nextflow.config` for test-specific settings:

```groovy
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

params {
    // Base directory for nf-core/modules test data
    modules_testdata_base_path   = 's3://ngi-igenomes/testdata/nf-core/modules/'
    // Base directory for pipeline-specific test data
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/PIPELINE_NAME/'

    // Input data
    input       = "${projectDir}/assets/samplesheet.csv"
    fasta       = "${params.pipelines_testdata_base_path}/reference/genome.fa.gz"
    fasta_index = "${params.pipelines_testdata_base_path}/reference/genome.fa.fai"
    outdir      = 'results'
}

// Impose sensible resource limits for testing
process {
    resourceLimits = [
        cpus: 2,
        memory: '3.GB',
        time: '2.h'
    ]
}

aws.client.anonymous = true // fixes S3 access issues on self-hosted runners

// Impose same minimum Nextflow version as the pipeline for testing
manifest {
    nextflowVersion = '!>=24.10.2'
}

// Disable all Nextflow reporting options
timeline { enabled = false }
report   { enabled = false }
trace    { enabled = false }
dag      { enabled = false }
``` 