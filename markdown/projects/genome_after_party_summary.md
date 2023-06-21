# sanger-tol/insdcdownload {#insdcdownload}

This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, lepidoptera, and requests) assemblies.

### Current features:
- Download genomes from NCBI.
- Unmask genome assembly.
- Build `samtools` `faidx` and `dict` indices for genome assembly.
- Create BED file with the coordinates of the masked region.
- Build `tabix` index for BED file.
- Custom download location.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/insdcdownload). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3). 

# sanger-tol/ensemblrepeatdownload {#ensemblrepeatdownload}

This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, lepidoptera, and requests) assemblies.

### Current features:
- Download the masked FASTA file from Ensembl.
- Extract the coordinates of the masked regions into a BED file.
- Compress and index the BED file with bgzip and tabix.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/ensemblrepeatdownload). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# sanger-tol/ensemblgenedownload {#ensemblgenedownload}

This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, lepidoptera, and requests) assemblies.

### Current features:
- Download from Ensembl gene annotation in GFF3 format.
- Download from Ensembl gene sequences in FASTA format.
- Compress and index all sequences files with `bgzip`, `samtools faidx`, and `samtools dict`.
- Compress and index the annotation files with `bgzip` and `tabix`.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/ensemblrepeatdownload). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# sanger-tol/sequencecomposition {#sequencecomposition}

This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, lepidoptera, and requests) assemblies.

### Current features:
- Run `fasta_windows` on the genome FASTA file.
- Extract single-statistics bedGraph files from the multi-statistics outputs.
- Compress and index all bedGraph and TSV files with `bgzip` and `tabix`.

### Planned features:
- Add [TRASH (Tandem Repeat Annotation and Structural Hierarchy)](https://github.com/vlothec/TRASH) pipeline
- Import features from [GDA (Genome Decomposition Analysis)](https://github.com/sanger-tol/gda) pipeline
- Update pipeline template
- Update samplesheet validation steps

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/sequencecomposition). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# sanger-tol/readmapping {#readmapping}

This pipeline is run for all **primary and alternative** ToL assemblies.

### Current features:
- Align short read data (HiC and Illumina) against the genome with `bwamem2 mem`
- Mark duplicates for short read alignment with `samtools`.
- Filter PacBio raw read data using vector database.
- Align long read data (ONT, PacBio CCS and PacBio CLR) against the genome with `minimap align`.
- Merge all alignment files at the individual level and convert to CRAM format.
- Calculate statistics for all alignment files using `samtools` `stats`, `flagstat`, and `idxstats`.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps
- Add calculation for how 
- Add support for Pacbio ULI reads
- Add support for RNAseq data
- Implement chunking to speed up alignment for all technologies
- Include additional metadata in aligned file headers
- Support compression with `crumble` for aligned files
- Support multiple output options – BAM, compressed BAM, CRAM, compressed CRAM

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/readmapping). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# sanger-tol/variantcalling {#variantcalling}

This pipeline is run for all **primary** ToL assemblies.

### Current features:
- Calls variants using DeepVariant for PacBio long read data.
- Speed improvements made by splitting the genome before calling variants.
- Outputs both VCF and GVCF formats.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps
- Add calculation for hetrozygosity, flag for homozygous alternatives, and InDel size distribution
- Create bedGraph for distribution of heterozygous sites across genome
- Add support for creating indices for aligned reads and the genome
- Add support for combining multiple libraries from the same sample
- Add support to map PacBio reads before calling variants
- Add structural variation detection
- Add variant calling for short read data with FreeBayes

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/variantcalling). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3). 

# sanger-tol/blobtoolkit {#blobtoolkit}

This pipeline will be run for all **primary** ToL assemblies after 1.0.0 release. Currently, the [Snakemake pipeline](https://github.com/blobtoolkit/blobtoolkit/tree/main/src/blobtoolkit-pipeline/src) is used in production.

### Current features:
- Calculate sequence statistics in 1kb windows for each contig.
- Count BUSCOs in 1kb windows for each contig using specific and basal lineages.
- Calculate coverage in 1kb windows using `blobtk depth`.
- Aggregate 1kb values into windows of fixed proportion (10%, 1% of contig length) and fixed length (100kb, 1Mb).
- `Diamond blastp` search of BUSCO gene models for basal lineages (archaea_odb10, bacteria_odb10 and eukaryota_odb10) against the UniProt reference proteomes.
- Import analysis results into a BlobDir dataset.
- BlobDir validation and static image generation.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps
- Improved BlobDir validation
- `Diamond blastx` search of assembly contigs against the UniProt reference proteomes
- NCBI `blastn` search of assembly contigs with no `Diamond blastx` match against the NCBI nt database

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/blobtoolkit). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# sanger-tol/genomenote {#genomenote}

This pipeline is run for all **primary** ToL assemblies.

### Current features:
- Create HiC contact map and chromosomal grid using `Cooler`.
- Create summary table using (1) assembly information, statistics and chromosome details from NCBI `datasets`, (2) genome completeness from `BUSCO`, (3) consensus quality and k-mer completeness from `MerquryFK`, and (4) HiC primary mapped percentage from `samtools flagstat`.

### Planned features:
- Update pipeline template
- Update samplesheet validation steps
- Improve genome metadata fetching and processing
- Combine results and metadata with template XML for submission to F1000

For more details, visit the [pipeline website](https://pipelines.tol.sanger.ac.uk/genomenote). If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).
