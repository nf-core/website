# Introduction

Genome After Party is a suite of pipelines to standardise the downstream analyses performed on all genomes produced by the Tree of Life. These include:

- **[sanger-tol/insdcdownload](/insdcdownload)** downloads assemblies from INSDC into a Tree of Life directory structure.
- **[sanger-tol/ensemblrepeatdownload](/ensemblrepeatdownload)** downloads repeat annotations from Ensembl into a Tree of Life directory structure.
- **[sanger-tol/ensemblgenedownload](/ensemblgenedownload)** downloads gene annotations from Ensembl into the Tree of Life directory structure.
- **[sanger-tol/sequencecomposition](/sequencecomposition)** extracts statistics from a genome about its sequence composition.
- **[sanger-tol/readmapping](/readmapping)** for aligning reads generated using Illumina, HiC, PacBio and Nanopore technologies against a genome assembly.
- **[sanger-tol/variantcalling](/variantcalling)** for calling variants using DeepVariant with PacBio data.
- **[sanger-tol/blobtoolkit](/blobtoolkit)** is used to identify and analyse non-target DNA for eukaryotic genomes.
- **[sanger-tol/genomenote](/genomenote)** creates HiC contact maps and collates (1) assembly information, statistics and chromosome details, (2) PacBio consensus quality and k-mer completeness, and (3) HiC mapping statistics.

Learn more about our pipelines on their dedicated pages. These pipelines are created using [Nextflow DSL2](https://www.nextflow.io) and [nf-core](https://nf-co.re) template. They are designed for portability, scalability and biodiversity.

A portal is being developed to automate the production of genome note publications. It will execute the Nextflow pipeline and populate an associated database with generated statistics and images. The platform is being designed in collaboration with the Enabling Platforms team to create genome note style publications for both internal Tree of Life assemblies as well as external genome assemblies.

> You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).
> Currently, we plan to **download** all primary and alternative ToL assemblies, as well as non-ToL (VGP, _Lepidoptera_, and requests) assemblies, and run the analysis pipelines **only on the primary** assemblies.
> [Let us know](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+) if some analyses would be useful to have on the other assemblies too.

# INSDC Download {#insdcdownload}

**[sanger-tol/insdcdownload](/insdcdownload)** downloads assemblies from INSDC into a Tree of Life directory structure. This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, _Lepidoptera_, and requests) assemblies.

### Current features:

- Download genome from NCBI as Fasta.
- Put the unmasked version under `assembly/release/` and the masked version under `analysis/`.
- Build `samtools` `faidx` and `dict` indices on the genome assemblies.
- Create BED file with the coordinates of the masked region.
- Compress and index the BED file with `bgzip` and `tabix`.

### Planned features:

- Prepare a file that can be used to populate BAM headers
- Update pipeline template
- Update samplesheet validation steps

If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Ensembl Repeat Download {#ensemblrepeatdownload}

**[sanger-tol/ensemblrepeatdownload](/ensemblrepeatdownload)** downloads repeat annotations from Ensembl into a Tree of Life directory structure. This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, _Lepidoptera_, and requests) assemblies.

### Current features:

- Download the masked FASTA file from Ensembl.
- Extract the coordinates of the masked regions into a BED file.
- Compress and index the BED file with `bgzip` and `tabix`.

### Planned features:

- Repeat density sub-workflow
- Update pipeline template
- Update samplesheet validation steps

If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Ensembl Gene Download {#ensemblgenedownload}

**[sanger-tol/ensemblgenedownload](/ensemblgenedownload)** downloads gene annotations from Ensembl into the Tree of Life directory structure. This pipeline is run for all **primary and alternative** ToL assemblies, as well as non-ToL (VGP, _Lepidoptera_, and requests) assemblies.

### Current features:

- Download from Ensembl gene annotation in GFF3 format.
- Download from Ensembl gene sequences in FASTA format.
- Compress and index all sequences files with `bgzip`, `samtools faidx`, and `samtools dict`.
- Compress and index the annotation files with `bgzip` and `tabix`.

### Planned features:

- Gene density sub-workflow
- Update pipeline template
- Update samplesheet validation steps

If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Sequence Composition {#sequencecomposition}

**[sanger-tol/sequencecomposition](/sequencecomposition)** extracts statistics from a genome about its sequence composition. This pipeline is run for all **primary** ToL assemblies.

### Current features:

- Run `fasta_windows` on the genome FASTA file.
- Extract single-statistics `bedGraph` files from the multi-statistics outputs.
- Compress and index all `bedGraph` and TSV files with `bgzip` and `tabix`.

### Planned features:

- Add [TRASH (Tandem Repeat Annotation and Structural Hierarchy)](https://github.com/vlothec/TRASH) pipeline
- Add [Pantera](https://github.com/piosierra/pantera) pipeline
- Add [EarlGreyTE](https://github.com/TobyBaril/EarlGrey) pipeline
- Add [stainedglass](https://github.com/mrvollger/StainedGlass) and/or [ModDotPlot](https://github.com/marbl/ModDotPlot) pipeline
- Convert all outputs to bigBed and build a public track-hub
- Update pipeline template
- Update samplesheet validation steps


### Potential features:

_Only if there is an actual demand! Use the form below to indicate interest._

- Import features from the [GDA (Genome Decomposition Analysis)](https://github.com/sanger-tol/gda) pipeline, for instance:
  - Low complexity repeats from `Dustmasker`
  - Inverted repeats from `einverted`
  - LTR retrotransposons from `LTRharvest` and `LTRdigest`
  - Tandem repeats from `trf`
  - Mappability track
  - Telomeric repeat annotation (tool to be confirmed)
  - Centromeric repeat annotation (tool to be confirmed)

If you have an idea for a new feature or would like this pipeline to run on other assemblies – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Read Mapping {#readmapping}

**[sanger-tol/readmapping](/readmapping)** aligns reads generated using Illumina, HiC, PacBio and Nanopore technologies against a genome assembly. This pipeline is run for all **primary** ToL assemblies.

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
- Add calculation for PacBio filtered data percentage
- Add support for Pacbio ULI reads
- Add support for RNAseq data
- Implement chunking to speed up alignment for all technologies
- Include additional metadata in aligned file headers
- Support compression with `crumble` for aligned files
- Support multiple output options – BAM, compressed BAM, CRAM, compressed CRAM

If you have an idea for a new feature or would like this pipeline to run on other assemblies – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Variant Calling {#variantcalling}

**[sanger-tol/variantcalling](/variantcalling)** calls (short) variants on PacBio data using DeepVariant. This pipeline is run for all **primary** ToL assemblies.

### Current features:

- Can combine multiple libraries from the same sample
- Optional read mapping subworkflow
- Calls variants using DeepVariant for PacBio long read data.
- Speed improvements made by splitting the genome before calling variants.
- Outputs both VCF and GVCF formats.

### Planned features:

- Update pipeline template
- Update samplesheet validation steps
- Add calculation for heterozygosity, flag for homozygous alternatives, and InDel size distribution
- Create bedGraph for distribution of heterozygous sites across genome
- Add support for creating indices for aligned reads and the genome
- Add structural variation detection
- Add variant calling for short read data
- Ability to process and analyse population genomics data

If you have an idea for a new feature or would like this pipeline to run on other assemblies – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# BlobToolKit {#blobtoolkit}

**[sanger-tol/blobtoolkit](/blobtoolkit)** is used to identify and analyse non-target DNA for eukaryotic genomes. This pipeline will be run for all **primary** ToL assemblies after 1.0.0 release. Currently, the [Snakemake version](https://github.com/blobtoolkit/blobtoolkit/tree/main/src/blobtoolkit-pipeline/src) is used in production.

### Current features:

- Calculate sequence statistics in 1kb windows for each contig.
- Count BUSCOs in 1kb windows for each contig using specific and basal lineages.
- Calculate coverage in 1kb windows using `blobtk depth`.
- Aggregate 1kb values into windows of fixed proportion (10%, 1% of contig length) and fixed length (100kb, 1Mb).
- `Diamond blastp` search of BUSCO gene models for basal lineages (archaea_odb10, bacteria_odb10 and eukaryota_odb10) against the UniProt reference proteomes.
- `Diamond blastx` search of assembly contigs against the UniProt reference proteomes
- NCBI `blastn` search of assembly contigs with no `Diamond blastx` match against the NCBI nt database
- Optional read mapping subworkflow
- Import analysis results into a BlobDir dataset.
- BlobDir validation and static image generation.

### Planned features:

- Update pipeline template

If you have an idea for a new feature or would like this pipeline to run on other assemblies – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).

# Genome Note {#genomenote}

**[sanger-tol/genomenote](/genomenote)** generates all the data (tables and figures) used in genome note publications. These include (1) assembly information, statistics and chromosome details, (2) PacBio consensus quality and k-mer completeness, and (3) HiC contact maps and mapping statistics. This pipeline is run for all **primary** ToL assemblies.

### Current features:

- Create HiC contact map and chromosomal grid using `Cooler`.
- Create summary table using (1) assembly information, statistics and chromosome details from NCBI `datasets`, (2) genome completeness from `BUSCO`, (3) consensus quality and k-mer completeness from `MerquryFK`, and (4) HiC primary mapped percentage from `samtools flagstat`.

### Planned features:

- Update pipeline template
- Update samplesheet validation steps
- Improve genome metadata fetching and processing
- Combine results and metadata with template Word document
- Add optional read mapping subworkflow
- Support Merqury and Meryl for k-mer analysis
- Run GFA stats

If you have an idea for a new feature or would like this pipeline to run on other assemblies – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+). You can see all planned features and requests on the [project board](https://github.com/orgs/sanger-tol/projects/3).
