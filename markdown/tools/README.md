<img align="left" src="https://goat.genomehubs.org/static/092eb939ef98ce8ebf19/about/browse.png" height="50px" style="margin-right: 8px">

## Genomes on a Tree (GoaT)

Genomes on a Tree (GoaT), an Elasticsearch-powered datastore and search index for genome-relevant metadata and sequencing project plans and statuses. Attributes can be queried against assembly and taxon indexes through 1) a mature [API](https://goat.genomehubs.org/api-docs), 2) a [web](https://goat.genomehubs.org) front end, and 3) a command line interface, [GoaT-cli](https://github.com/genomehubs/goat-cli). GoaT uses NCBI Taxonomy as backbone and allows retrieval of inferred estimates for missing attribute values based on phylogenetic interpolation. The web front end additionally provides summary visualisations for data exploration and reporting (see https://goat.genomehubs.org). GoaT was implemented as an instance of the open-source GenomeHubs codebase. For more information visit https://github.com/genomehubs/genomehubs.

**Pipelines:** [genomenote](https://pipelines.tol.sanger.ac.uk/genomenote) | [blobtoolkit](https://pipelines.tol.sanger.ac.uk/blobtoolkit)


<img align="left" src="assets/img/Cobiont_sm.png" height="50px" style="margin-right: 8px">

## [CobiontID](https://cobiontid.github.io/) – Identifying cobionts and contaminants

The **CobiontID** process has two parts: First, **Marker scan** provides taxonomic information. HMM profiles of marker genes, such as rRNAs, which are well-sampled and conserved, are useful to classify sequences from genomes that are otherwise too diverged from their closest sequenced relative. We can therefore gauge which species are present in a given sample, and construct streamlined databases for read classification. Second, a combination of assembly, read mapping and **compositional clustering** allows the sequences to be assigned to groups that can be tagged with this taxonomic information.

**Pipelines:** cobiontcheck


## [YaHS](http://dx.doi.org/10.1093/bioinformatics/btac808) – yet another Hi-C scaffolding tool

YaHS is a scaffolding tool using Hi-C data. It relies on a new algorithm for contig joining detection which considers the topological distribution of Hi-C signals aiming to distinguish real interaction signals from mapping noises. YaHS has been tested in a wide range of genome assemblies. Compared to other Hi-C scaffolding tools, it usually generates more contiguous scaffolds - especially with a higher N90 and L90 statistics. It is also super fast - takes less than 5 minutes to reconstruct the human genome from an assembly of 5,483 contigs with ~45X Hi-C data.

**Pipelines:** [genomeassembly](https://pipelines.tol.sanger.ac.uk/genomeassembly)


## [MitoHiFi](https://www.biorxiv.org/content/10.1101/2022.12.23.521667v2) – a python pipeline for mitochondrial genome assembly from PacBio HiFi reads

MitoHiFi is able to assemble mitochondrial genomes from a wide phylogenetic range of taxa from Pacbio HiFi data. MitoHiFi is written in python and is freely available on GitHub. MitoHiFi has been used to assemble 374 mitochondrial genomes (369 from 12 phyla and 39 orders of Metazoa and from 6 species of Fungi) for the Darwin Tree of Life Project, the Vertebrate Genomes Project and the Aquatic Symbiosis Genome Project

**Pipelines:** [genomeassembly](https://pipelines.tol.sanger.ac.uk/genomeassembly)
