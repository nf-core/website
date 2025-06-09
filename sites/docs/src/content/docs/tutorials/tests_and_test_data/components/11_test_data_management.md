---
title: "11. Test Data Management"
subtitle: Organizing and managing test datasets
weight: 110
---

## Test Data Strategy

### nf-core Test Data Repository

nf-core maintains centralized test datasets:

- **Modules test data**: `https://github.com/nf-core/test-datasets/tree/modules`
- **Pipeline test data**: `https://github.com/nf-core/test-datasets/tree/<pipeline-name>`

```groovy
// Standard nf-core test data reference
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
```

### Local Test Data

For component-specific or sensitive data:

```
tests/
├── data/
│   ├── genomics/
│   │   ├── reference.fa
│   │   ├── sample1.fastq.gz
│   │   └── sample2.fastq.gz
│   └── proteomics/
│       ├── proteins.fasta
│       └── spectra.mzML
├── samplesheets/
│   ├── minimal.csv
│   └── full.csv
└── configs/
    ├── test.config
    └── test_full.config
```

## Data Organization Patterns

### By Data Type

```
test-datasets/
├── genomics/
│   ├── homo_sapiens/
│   │   ├── genome/
│   │   │   ├── genome.fasta
│   │   │   └── genome.fasta.fai
│   │   ├── illumina/
│   │   │   ├── bam/
│   │   │   ├── fastq/
│   │   │   └── vcf/
│   │   └── annotations/
│   │       ├── genes.gtf
│   │       └── repeats.bed
│   └── sarscov2/
│       ├── genome/
│       └── illumina/
└── proteomics/
    ├── human/
    └── ecoli/
```

### By Size/Complexity

```
test-data/
├── minimal/          # Ultra-small datasets for CI
│   ├── reads_10k.fastq.gz
│   └── ref_100bp.fa
├── small/           # Small representative datasets
│   ├── reads_100k.fastq.gz
│   └── ref_1kb.fa
├── medium/          # Realistic but manageable datasets
│   ├── reads_1M.fastq.gz
│   └── ref_100kb.fa
└── large/           # Full-scale test datasets
    ├── reads_10M.fastq.gz
    └── ref_1Mb.fa
```

## Test Data Access Patterns

### Remote Data Access

```groovy
// Using parameters for flexible data sources
nextflow_process {
    name "Test Process FASTQC"
    script "../main.nf"
    process "FASTQC"

    test("Remote test data") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test' ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
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

### Local Data with Fallbacks

```groovy
test("Local data with remote fallback") {
    when {
        process {
            """
            def localFile = file('tests/data/test_sample.fastq.gz')
            def remoteFile = file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
            def testFile = localFile.exists() ? localFile : remoteFile

            input[0] = [
                [ id:'test' ],
                testFile
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

### Dynamic Data Generation

```groovy
test("Generated test data") {
    setup {
        // Generate test data on-the-fly
        def testFastq = new File("tests/generated_test.fastq")
        testFastq.text = """@read1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@read2
GCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIII
"""
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/generated_test.fastq', checkIfExists: true)
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

## Data Size Management

### Minimal Test Data Creation

```bash
# Create minimal FASTQ files
head -n 400 large_sample.fastq > minimal_sample.fastq
gzip minimal_sample.fastq

# Create minimal reference
head -n 100 reference.fa > minimal_reference.fa

# Create minimal BAM files
samtools view -s 0.01 large.bam | samtools view -bS - > minimal.bam
```

### Data Subsetting Strategies

```groovy
// Function to create test data subsets
def createMinimalDataset(originalFile, outputFile, fraction = 0.001) {
    if (originalFile.name.endsWith('.fastq.gz')) {
        // Sample FASTQ reads
        """
        zcat ${originalFile} | head -n 4000 | gzip > ${outputFile}
        """.execute()
    } else if (originalFile.name.endsWith('.bam')) {
        // Sample BAM alignments
        """
        samtools view -s ${fraction} -b ${originalFile} > ${outputFile}
        """.execute()
    }
}

test("With minimal dataset") {
    setup {
        createMinimalDataset(
            file('full_dataset.fastq.gz'),
            file('tests/minimal_dataset.fastq.gz')
        )
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/minimal_dataset.fastq.gz', checkIfExists: true)
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

## Environment-Specific Data Management

### CI/CD Optimized Data

```groovy
profiles {
    ci {
        params {
            // Use minimal datasets for CI
            modules_testdata_base_path = 'https://github.com/nf-core/test-datasets/raw/modules/data/'
            max_memory = '6.GB'
            max_cpus = 2
        }
    }

    local_dev {
        params {
            // Use local mirror for development
            modules_testdata_base_path = '/data/test-datasets/modules/'
            max_memory = '32.GB'
            max_cpus = 8
        }
    }
}
```

### Data Caching Strategies

```groovy
// Cache test data locally
def cacheTestData(remoteUrl, localPath) {
    def localFile = file(localPath)
    if (!localFile.exists()) {
        localFile.parentFile.mkdirs()
        new URL(remoteUrl).withInputStream { input ->
            localFile.withOutputStream { output ->
                output << input
            }
        }
    }
    return localFile
}

test("Cached data access") {
    when {
        process {
            """
            def cachedFile = cacheTestData(
                '${params.modules_testdata_base_path}genomics/sarscov2/illumina/fastq/test_1.fastq.gz',
                'tests/cache/test_1.fastq.gz'
            )

            input[0] = [
                [ id:'test' ],
                cachedFile
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

## Data Validation and Quality

### Data Integrity Checks

```groovy
test("Data validation") {
    setup {
        // Validate test data before use
        def validateFastq = { file ->
            def lines = file.readLines()
            assert lines.size() % 4 == 0, "Invalid FASTQ format"
            assert lines[0].startsWith("@"), "Missing FASTQ header"
            assert lines[2].startsWith("+"), "Missing FASTQ separator"
            return true
        }

        def testFile = file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz')
        validateFastq(testFile)
    }

    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
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

### Format-Specific Validation

```groovy
// Custom validation functions
def validateBAM(bamFile) {
    def header = "samtools view -H ${bamFile}".execute().text
    assert header.contains("@HD"), "BAM missing header"
    assert header.contains("SO:"), "BAM missing sort order"
    return true
}

def validateFASTA(fastaFile) {
    def content = fastaFile.text
    assert content.contains(">"), "FASTA missing header"
    assert content.matches(/(?s)>[^\n]+\n[ATCGN\n]+/), "Invalid FASTA format"
    return true
}

def validateVCF(vcfFile) {
    def lines = vcfFile.readLines()
    assert lines.any { it.startsWith("##fileformat=VCF") }, "VCF missing format line"
    assert lines.any { it.startsWith("#CHROM") }, "VCF missing header line"
    return true
}
```

## Synthetic Data Generation

### Creating Realistic Test Data

```groovy
// Generate synthetic genomic data
def generateSyntheticReads(outputFile, numReads = 1000, readLength = 150) {
    def bases = ['A', 'T', 'C', 'G']
    def qualities = (20..40).collect { (char)(it + 33) }

    outputFile.withWriter { writer ->
        (1..numReads).each { i ->
            // FASTQ header
            writer.println("@synthetic_read_${i}")

            // Sequence
            def sequence = (1..readLength).collect {
                bases[new Random().nextInt(4)]
            }.join('')
            writer.println(sequence)

            // Separator
            writer.println("+")

            // Quality scores
            def quality = (1..readLength).collect {
                qualities[new Random().nextInt(qualities.size())]
            }.join('')
            writer.println(quality)
        }
    }
}

test("Synthetic data test") {
    setup {
        generateSyntheticReads(file('tests/synthetic_reads.fastq'))
    }

    when {
        process {
            """
            input[0] = [
                [ id:'synthetic' ],
                file('tests/synthetic_reads.fastq', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match("synthetic_data") }
        )
    }
}
```

## Data Privacy and Security

### Sensitive Data Handling

```groovy
// Configuration for sensitive data
profiles {
    sensitive_data {
        params {
            // Use encrypted or anonymized datasets
            test_data_base = '/secure/anonymized-datasets/'

            // Ensure outputs are properly handled
            publish_dir_mode = 'copy'
            cleanup_intermediate = true
        }

        process {
            // Add security labels
            label = 'sensitive_data'

            // Ensure cleanup
            cleanup = true
        }
    }
}
```

### Data Anonymization

```groovy
def anonymizeData(inputFile, outputFile) {
    // Remove or replace sensitive identifiers
    def content = inputFile.text
    content = content.replaceAll(/patient_\d+/, 'sample_ID')
    content = content.replaceAll(/\b\d{4}-\d{2}-\d{2}\b/, 'YYYY-MM-DD')
    outputFile.text = content
}

test("Anonymized data test") {
    setup {
        anonymizeData(
            file('sensitive/patient_data.txt'),
            file('tests/anonymized_data.txt')
        )
    }

    when {
        process {
            """
            input[0] = [
                [ id:'anonymized' ],
                file('tests/anonymized_data.txt', checkIfExists: true)
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

## Best Practices

### 1. Use Appropriate Data Sizes

- **CI tests**: Ultra-minimal datasets (< 1MB)
- **Development**: Small representative datasets (1-10MB)
- **Integration**: Realistic datasets (10-100MB)
- **Performance**: Full-scale datasets (> 100MB)

### 2. Version Control for Test Data

```bash
# Use Git LFS for larger test files
git lfs track "*.fastq.gz"
git lfs track "*.bam"
git lfs track "*.vcf.gz"

# Or use external storage with checksums
echo "test_data.fastq.gz  md5:1234567890abcdef" > checksums.md5
```

### 3. Document Data Sources

```groovy
// Document test data provenance
params {
    test_data_info = [
        'test_sample.fastq.gz': [
            source: 'SRA:SRR123456',
            description: 'E. coli genome sequencing, paired-end reads',
            size: '2.5MB',
            checksum: 'md5:abcdef1234567890'
        ]
    ]
}
```

## Next Steps

Continue to [Error Testing & Edge Cases](./12_error_testing_edge_cases.md) to learn about testing failure scenarios and edge cases.
