---
title: Data management
subtitle: How to plan your project, estimate resources, and share your results.
---

## Data management

Funding agencies are recognizing the importance of research data management and some now request detailed Data Management Plans (DMP) as part of the grant application. Research data management concerns the organization, storage, preservation, and sharing of data that is collected or analyzed during a research project. Proper planning and data management facilitates sharing and allows others to validate and reuse the data. Guidance is provided below to aid the creation of DMPs, estimate resources needed by nf-core workflows, and how to share the resulting data.

### Data Management Plan

A Data Management Plan (DMP) is a revisable document explaining how you intend to handle new and existing data, during and following the conclusion of your research project.
It is wise to write a DMP as early as possible, using either a tool provided by your host institution or for example
[DS Wizard](https://ds-wizard.org/) or [DMP Online](https://dmponline.dcc.ac.uk/).
Ethical and legal considerations regarding the data will depend on where the research is conducted, this is especially true for projects including sensitive human data. For more information about the Swedish context, please review this page on [Sensitive personal data](https://scilifelab-data-guidelines.readthedocs.io/en/latest/docs/general/sensitive_data.html).

### Data storage and computational resources

To estimate computational resources needed for a specific pipeline please see the pipeline documentation. This lists the different output filetypes you can expect.
In the future we hope to automate a full run of each pipeline after every release. The pipeline docs will then show a full set of results from a real run, along with all file sizes. This can then be used as a guide as to what to expect for your data.

Backing up and archiving your data is essential. The 3-2-1 rule of thumb means that you should have 3 copies of the data, on 2 different types of media, and 1 of the copies at different physical location. Consider uploading the raw data to a repository already when receiving them, under an embargo (if that is important to you). This way you always have an off-site backup with the added benefit of making the data sharing phase more efficient. Identifying a suitable repository early on will allow you to conform to their standards and metadata requirements already from the start.

Archiving is often the responsibility of your host institution, contact them for more details.

### Data compression

Some tools require compressed input files, which have many advantages: they take less space for storage and sharing. The most frequent format is gzip; it is accepted by many tools which you can check in tool manuals. To compress a file, you can use [bzip2](https://sourceware.org/bzip2/) which creates a non-blocked compressed file. If a tool only accepts uncompressed file input, you can uncompress the file and parse it via a pipe to the tool without saving the compressed version of the input. Here is an example:

```
gzip input | TOOL > output

```

If a tool requires a blocked compressed file (BGZP), in which the information is more easily accessible than a non-blocked compression format, you can use the [htslib/bgzip](http://www.htslib.org/doc/bgzip.html) tool. This is typically needed by SAMTOOLS during sequence alignment analyses.

Read a more detailed explanation on which format to choose [here](https://www.uppmax.uu.se/support/faq/resources-faq/which-compression-format-should-i-use-for-ngs-related-files/).

## Data sharing

In the era of [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable and Reusable) and [Open science](https://ec.europa.eu/research/openscience/index.cfm), datasets should be made available to the public, for example by submitting your data to a public repository.

### Choosing a repository

It’s recommended to choose a domain-specific repository when possible. It is also important to consider the sustainability of the repository to ensure that the data will remain public. Please see the [EBI archive wizard](https://www.ebi.ac.uk/submission/) or [SciLifeLab's data guidelines](https://scilifelab-data-guidelines.readthedocs.io/en/latest/docs/index.html) for suggestions depending on data type. You can also refer to the [ELIXIR Deposition Databases](https://elixir-europe.org/services/tag/elixir-deposition-databases), [Scientific Data’s Recommended Data Repositories](https://www.nature.com/sdata/policies/repositories), [FAIRsharing.org](https://fairsharing.org/databases/) and [re3data.org](https://www.re3data.org/) to find suitable repositories.
Also note that funding agencies might have specific requirements regarding data deposition. For example, data generated in projects funded by US federal grants should be deposited into public databases such as [SRA](https://www.ncbi.nlm.nih.gov/sra) for raw sequencing data and [GEO](https://www.ncbi.nlm.nih.gov/geo/) for functional genomics data.

For datasets that do not fit into domain-specific repositories, you can use an institutional repository when available or a general repository such as [Figshare](https://figshare.com/) and [Zenodo](https://zenodo.org/).

### Preparing for submission

#### Describing and organizing your data

Metadata should be provided to help others discover, identify and interpret the data. Researchers are strongly encouraged to use community metadata standards and ontologies where these are in place, consult e.g [FAIRsharing.org](https://fairsharing.org/standards/). Data repositories may also provide guidance about metadata standards and requirements. Capture any additional documentation needed to enable reuse of the data in Readme text files and [Data Dictionaries](https://help.osf.io/hc/en-us/articles/360019739054-How-to-Make-a-Data-Dictionary) that describe what all the variable names and values in your data really mean. Identifiers to refer to e.g. ontology terms can be designed for computers or for people; in a FAIR data context it is recommended to supply both a human-readable as well as a machine-resolvable Persistent Identifier (PID) for each concept used in the data.

#### Data integrity

Some repositories require md5 checksums to be uploaded along with the files. Also consider adding checks that your data files follow the intended file formats and can be opened by standard software for those formats.

If you are using a Linux system, you can generate md5 checksums using the `md5sum` command.

#### Choosing a license

To ascertain re-usability data should be released with a clear and accessible data usage license. We suggest making your data available under licenses that permit free reuse of data, e.g. a Creative Commons license, such as CC0 or CC-BY. The [EUDAT license selector wizard](https://ufal.github.io/public-license-selector/) can help you select suitable licenses for your data.
Note that sequence data submitted to [ENA](https://www.ebi.ac.uk/ena)/[GenBank](https://www.ncbi.nlm.nih.gov/genbank/)/[DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) are implicitly free to reuse by others as specified in the [INCD Standards and policies](https://www.ebi.ac.uk/ena/standards-and-policies).
