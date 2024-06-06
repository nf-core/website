Do you know of a ToL tool that we are missing? [Let us know.](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato&labels=tool%2Cenhancement&projects=&template=add_tool.yaml&title=%5BTool%5D%3A+)

<img align="left" src="/assets/img/genomehubs-icon-512.png" height="50px" style="margin-right: 8px">

## Genomes on a Tree (GoaT)

Genomes on a Tree (GoaT), an Elasticsearch-powered datastore and search index for genome-relevant metadata and sequencing project plans and statuses. Attributes can be queried against assembly and taxon indexes through 1) an [API](https://goat.genomehubs.org/api-docs), 2) a [web](https://goat.genomehubs.org) front end, and 3) a command line interface, [GoaT-cli](https://github.com/genomehubs/goat-cli). GoaT uses NCBI Taxonomy as backbone and allows retrieval of inferred values for missing attribute values based on phylogenetic interpolation. The web front end additionally provides summary visualisations for data exploration and reporting (see https://goat.genomehubs.org). GoaT was implemented as an instance of the open-source GenomeHubs codebase. For more information visit https://github.com/genomehubs/genomehubs.

**Pipelines:** [genomenote](/genomenote) | [blobtoolkit](/blobtoolkit)

</br>

<img align="left" src="/assets/img/cobiontID.png" height="50px" style="margin-right: 8px">

## CobiontID

[CobiontID](https://cobiontid.github.io/) identifies cobionts and contaminants in two stages. First, marker scan provides taxonomic information. HMM profiles of marker genes, such as rRNAs, which are well-sampled and conserved, are useful to classify sequences from genomes that are otherwise too diverged from their closest sequenced relative. We can therefore gauge which species are present in a given sample, and construct streamlined databases for read classification. Second, a combination of assembly, read mapping and compositional clustering allows the sequences to be assigned to groups that can be tagged with this taxonomic information.

**Pipelines:** cobiontcheck

</br>

## Yet another Hi-C Scaffolding tool (YaHS)

[YaHS](http://dx.doi.org/10.1093/bioinformatics/btac808) is a scaffolding tool using Hi-C data. It relies on a new algorithm for contig joining detection which considers the topological distribution of Hi-C signals aiming to distinguish real interaction signals from mapping noises. YaHS has been tested in a wide range of genome assemblies. Compared to other Hi-C scaffolding tools, it usually generates more contiguous scaffolds - especially with a higher N90 and L90 statistics. It is also super fast - takes less than 5 minutes to reconstruct the human genome from an assembly of 5,483 contigs with ~45X Hi-C data.

**Pipelines:** [genomeassembly](/genomeassembly)

</br>

## MitoHiFi

[MitoHiFi](https://www.biorxiv.org/content/10.1101/2022.12.23.521667v2) is able to assemble mitochondrial genomes from a wide phylogenetic range of taxa from Pacbio HiFi data. MitoHiFi is written in python and is freely available on GitHub. MitoHiFi has been used to assemble 374 mitochondrial genomes (369 from 12 phyla and 39 orders of Metazoa and from 6 species of Fungi) for the Darwin Tree of Life Project, the Vertebrate Genomes Project and the Aquatic Symbiosis Genome Project

**Pipelines:** [genomeassembly](/genomeassembly)

<br/>

## Oatk

Developed for _de novo_ assembly of complex plant organelle genomes using PacBio HiFi data, [Oatk](https://github.com/c-zhou/oatk) can also be used to assemble other, simpler organelle genomes such as animal mitochondria. From whole genome sequencing data input, Oatk's tools produce annotation (using profile HMMs from [OatkDB](https://github.com/c-zhou/OatkDB.git)) and cirularised graph assemblies of plastid and mitochondrial genomes.

**Pipelines:** [genomeassembly](/genomeassembly)

<br/>

## Pretext

Pretext is a suite of tools used to curate the chromosome level structure of long read assemblies of novel genomes using Hi-C contact maps:

- [PretextMap](https://github.com/sanger-tol/PretextMap) - Produces a Pretext map file from Hi-C data mapped to an assembly in SAM/BAM/CRAM or [Pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) format. The Pretext map file stores the Hi-C contact map in an OpenGL texture.

- [PretextGraph](https://github.com/sanger-tol/PretextGraph) - Used to add BED-Graph format anntoation to a Pretext map file.

- [PretextView](https://github.com/sanger-tol/PretextView) - An interactive, OpenGL based viewer and editor of Pretext map files, used to curate the chromosome level structure of assemblies. The edited genome is exported in an AGP format file, and includes tags introduced during curation which label contaminants, haplotigs, chromosome names and haplotypes.

- [PretextSnapshot](https://github.com/sanger-tol/PretextSnapshot) - Extracts Hi-C contact map images from Pretext map files.

**Pipelines:** [genomeassembly](/genomeassembly) | [treeval](/treeval)

<br/>

## AGP / TPF Utilities

Nucleotide coodinates (other than the start or end of contigs) in AGP files produced by [PretextView](https://github.com/sanger-tol/PretextView) are only accurate to the number of bases spanned by each pixel in its map file. The `pretext-to-tpf` command of [AGP / TPF Utilites](https://github.com/sanger-tol/agp-tpf-utils/) uses the orginal, unedited assembly to correct the coordinates of the assembly contigs. It also uses tags from PretextView to label chromosomes, and to separate haplotypes, haplotigs and contaminants.
