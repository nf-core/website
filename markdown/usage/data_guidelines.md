---
title: Data Management Guidelines
subtitle:
---

## Data management

### Data Management Plan

A Data Management Plan (DMP) is a revisable document explaining how you intend to handle new and existing data, during and following the conclusion of your research project. It is wise to write a DMP as early as possible, using either a tool provided by your host institution or [DS Wizard](http://dsw.scilifelab.se/) ([users outside Sweden](https://ds-wizard.org/)). Ethical and legal considerations regarding the data will depend on where the research is conducted, this is especially true for projects including sensitive human data. For more  information about the Swedish context, please review this page on [Sensitive personal data](https://scilifelab-data-guidelines.readthedocs.io/en/latest/docs/general/sensitive_data.html).

### Data storage and computational resources

To estimate computational resources needed for a specific pipeline see docs/output page

/// Insert statement/instructions on where or how to get the estimated file sizes, folder structures etc. ///

Backing up and archiving your data is essential. The 3-2-1 rule of thumb means that you should have 3 copies of the data, on 2 different types of media, and 1 of the copies at different physical location. Consider uploading the raw data to a repository already when receiving them, under an embargo (if that is important to you). This way you always have an off-site backup with the added benefit of making the data sharing phase more efficient. Identifying a suitable repository early on will allow you to conform to their standards and metadata requirements already from the start.

Archiving is often the responsibility of your host institution, contact them for more details.

## Data sharing

In the era of [FAIR](https://www.force11.org/group/fairgroup/fairprinciples) (Findable, Accessible, Interoperable and Reusable) and [Open science](https://www.vr.se/english/mandates/open-science/open-access-to-research-data.html), datasets should be made available to the public, for example by submitting your data to a public repository.

### Choosing a repository

It’s recommended to choose a domain-specific repository when possible. It is also important to consider the sustainability of the repository to ensure that the data will remain public. Please see[ SciLifeLab's data guidelines](https://scilifelab-data-guidelines.readthedocs.io/en/latest/docs/index.html) or the [EBI archive wizard](https://www.ebi.ac.uk/submission/) for suggestions depending on data type. You can also refer to the [ELIXIR Deposition Databases](https://elixir-europe.org/services/tag/elixir-deposition-databases) and [Scientific Data’s Recommended Data Repositories](https://www.nature.com/sdata/policies/repositories), to find suitable repositories.

For datasets that do not fit into domain-specific repositories, you can use an institutional repository when available (e.g. [SciLifeLab Data Repository](https://scilifelab.figshare.com/)) or a general repository such as [Figshare](https://figshare.com/) and [Zenodo](https://zenodo.org/).

### Preparing for submission

#### **Describing and organizing your data**

Metadata should be provided to help others discover, identify and interpret the data. Researchers are strongly encouraged to use community metadata standards and ontologies where these are in place, consult e.g [FAIRsharing.org](https://fairsharing.org/databases/). Data repositories may also provide guidance about metadata standards and requirements. Capture any additional documentation needed to enable reuse of the data in Readme text files and [Data Dictionaries](https://help.osf.io/hc/en-us/articles/360019739054-How-to-Make-a-Data-Dictionary) that describe what all the variable names and values in your data really mean. Identifiers to refer to e.g. ontology terms can be designed for computers or for people; in a FAIR data context it is recommended to supply both a human-readable as well as a machine-resolvable Persistent Identifier (PID) for each concept used in the data.

#### **Data integrity**

Some repositories require md5 checksums to be uploaded along with the files. Also consider adding checks that your data files follow the intended file formats and can be opened by standard software for those formats.

///Consult docs/output page to see if the pipeline generates the sums automatically and where to find them. ///

#### **Choosing a license**

To ascertain re-usability data should be released with a clear and accessible data usage license. We suggest making your data available under licences that permit free reuse of data, e.g. a Creative Commons licence, such as CC0 or CC-BY. The [EUDAT licence selector wizard ](https://ufal.github.io/public-license-selector/)can help you select suitable licences for your data. Note that sequence data submitted to ENA (or GenBank) are implicitly free to reuse by others as specified in the [INCD Standards and policies]( https://www.ebi.ac.uk/ena/standards-and-policies).
