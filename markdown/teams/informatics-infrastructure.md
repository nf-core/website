## Introduction {#introduction}

The [Tree of Life](https://www.sanger.ac.uk/programme/tree-of-life/) projects will generate tens of thousands of high-quality genomes – more than have ever been sequenced! It is a challenging and extremely exciting task that will shape the future of biology, and the team’s role is to provide the platform for assembling and analysing those genomes at an unprecedented scale. We are the interface between the Tree of Life teams (assembly production and faculty research) and Sanger IT, working together with the informatics teams of the other programmes.

The team is organised in three poles.

:open_file_folder: **Data management**: Our data curators and managers maintain the integrity, consistency, and quality, or multiple databases used in production, including Genomes on a Tree (GoaT), Sample Tracking System (STS), Collaborative Open Plant Omics (COPO), and BioSamples.

:octocat: **Bioinformatics**: Our bioinformaticians develop the suite of analysis pipelines that will run on every genome produced in Tree of Life, providing a central database of core results available for all.

:nut_and_bolt: **Systems**: We develop and maintain some core systems used in production, including the execution and tracking of all bioinformatics pipelines, and the deployment of third-party web applications for internal use.

## Tech stack {#tech-stack}

The team uses a wide range of technologies, frameworks and programming languages, including Nextflow, Python, Conda, Jira, LSF, Singularity, and Kubernetes. The technology wheel below shows most of their logos. How many can you recognise?

<img src="https://www.sanger.ac.uk/wp-content/uploads/Informatics-Infrastructure-Technology-wheel-300.png" alt="Tech stack" height=500px>

## Projects {#projects}

### Genome After Party {#genome-after-party}

🧮 **Tech Stack:** [Nextflow DSL2](https://www.nextflow.io), [Python](https://www.python.org), [React](https://react.dev), [SQLAlchemy](https://www.sqlalchemy.org), and [PostgreSQL](https://www.postgresql.org)

[Genome After Party](https://pipelines.tol.sanger.ac.uk/genome_after_party) is a suite of pipelines to standardise the downstream analyses performed on all genomes produced by the Tree of Life. These include:

- **[sanger-tol/insdcdownload](https://pipelines.tol.sanger.ac.uk/insdcdownload)** downloads assemblies from INSDC into a Tree of Life directory structure.
- **[sanger-tol/ensemblrepeatdownload](https://pipelines.tol.sanger.ac.uk/ensemblrepeatdownload)** downloads repeat annotations from Ensembl into a Tree of Life directory structure.
- **[sanger-tol/ensemblgenedownload](https://pipelines.tol.sanger.ac.uk/ensemblgenedownload)** downloads gene annotations from Ensembl into the Tree of Life directory structure.
- **[sanger-tol/sequencecomposition](https://pipelines.tol.sanger.ac.uk/sequencecomposition)** extracts statistics from a genome about its sequence composition.
- **[sanger-tol/readmapping](https://pipelines.tol.sanger.ac.uk/readmapping)** for aligning reads generated using Illumina, HiC, PacBio and Nanopore technologies against a genome assembly.
- **[sanger-tol/variantcalling](https://pipelines.tol.sanger.ac.uk/variantcalling)** for calling variants using DeepVariant with PacBio data.
- **[sanger-tol/blobtoolkit](https://pipelines.tol.sanger.ac.uk/blobtoolkit)** is used to identify and analyse non-target DNA for eukaryotic genomes.
- **[sanger-tol/genomenote](https://pipelines.tol.sanger.ac.uk/genomenote)** creates HiC contact maps and collates (1) assembly information, statistics and chromosome details, (2) PacBio consensus quality and k-mer completeness, and (3) HiC mapping statistics.

Learn more about our pipelines on their dedicated pages. These pipelines are created using [Nextflow DSL2](https://www.nextflow.io) and [nf-core](https://nf-co.re) template. They are designed for portability, scalability and biodiversity. 

A portal is being developed to automate the production of genome note publications. It will execute the Nextflow pipeline and populate an associated database with generated statistics and images. The platform is being designed in collaboration with the Enabling Platforms team to create genome note style publications for both internal Tree of Life assemblies as well as external genome assemblies.

If you have an idea for a new feature – [send us your request](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=priyanka-surana&labels=pipeline%2Cenhancement&projects=&template=genome_after_party_feature_request.yaml&title=%5BFeature%5D%3A+).

## Members {#members}

### Matthieu Muffato, Team Lead {#matthieu-muffato}

[<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/muffato-matthieu/) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/github.svg" height=25px>](https://github.com/muffato) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/matthieu-muffato/) &nbsp;

</br> Matthieu leads the Informatics Infrastructure team, which guides the implementation and delivery of the genome assembly pipelines, and provides support for large-scale genome analyses for the Tree of Life faculty teams. He joined the Wellcome Sanger Institute in February 2021, to form the Informatics Infrastructure team for the Tree of Life programme. He has recruited 7 team members, with skills covering data curation & management, software development & operations, and bioinformatics.

### Priyanka Surana, Senior Bioinformatician {#priyanka-surana}

[<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/surana-priyanka/) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/github.svg" height=25px>](https://github.com/priyanka-surana) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/priyanka-surana/) &nbsp;

</br> Priyanka is a Senior Bioinformatician, overseeing the development of Nextflow pipelines for genome assembly, curation and downstream analyses. She also facilitates the workflows community and is passionate about building networks that support peer learning.

### Cibele G Sotero-Caio, Genomic Data Curator {#cibele-sotero-caio}

[<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/sotero-caio-cibele) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/github.svg" height=25px>](https://github.com/ccaio) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/cibele-sotero-caio-b379071a6/) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/twitter.svg" height=25px>](https://twitter.com/CibeleCaio) &nbsp;

</br> Cibele is the data curator for the Genomes on a Tree (GoaT) - a platform developed to support the Tree of Life and other sequencing initiatives of the Earth Biogenome project (EBP).

### Paul Davis, Data Manager {#paul-davis}

[<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/davis-paul/) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/github.svg" height=25px>](https://github.com/Paul-Davis) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/paul-davis-uk/) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/twitter.svg" height=25px>](https://twitter.com/SirPaulDavis) &nbsp;

</br> Paul works on the main ToL Genome Engine. This system was developed by the ToL to manage and track samples from collection, onboarding, processing in the lab, sequencing and finally the publication of the assembly and Genome Note publication. As there are many steps in this process developing methodology to identify issues as early as possible is vital to avoid wasted time and resource. Paul works at all levels of the project fielding questions about data flow, data fixes and helps other ToL staff and project stakeholders with data and information. Paul also interacts with external groups and stakeholders to maintain data integrity in the public domain.

### Beth Yates, Bioinformatics Engineer {#beth-yates}

[<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/github.svg" height=25px>](https://github.com/BethYates) &nbsp; [<img align="left" src="https://raw.githubusercontent.com/sanger-tol/pipelines-website/main/public_html/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/bethanyates/) &nbsp;

</br> Beth is a Bioinformatics Engineer working on a building a platform to automate the production of [Genome Note publications](https://wellcomeopenresearch.org/treeoflife). The Universal Genome Note platform consists of a web portal, database and Nextflow pipelines. Beth is contributing to the genomenote pipeline, this pipeline fetches assembly meta data and generates some of the figures and statistics included in each genome note.

### Alumni {#alumni}

- Zaynab Butt, Informatics and Digital Associate
