## Introduction

The sanger-tol project is the informatics branch of the Tree of Life programme at the Wellcome Sanger Institute, UK.
Here you can read how we organise ourselves, how we are funded and how the sanger-tol project was started.

We welcome external contributions and collaboration on our projects.
If you'd like to be involved, [drop us a message](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato%2Cmuffato&labels=connect&projects=&template=contact_us.yaml&title=%5BContact+Us%5D%3A+) or comment directly on an existing GitHub issue.

Please note that all community members are expected to adhere to our [code of conduct](/code_of_conduct), which is adopted from [nf-core](https://nf-co.re).

## Informatics Infrastructure team {#it}

### Introduction {#it-intro}

The Informatics Infrastructure team provides support for the production of reference genome assemblies and large-scale genome analyses in the Tree of Life programme, and helps with the management and use of informatics resources and digital solutions.

The team is organised in three poles.

📂 **Data management**: Our data curators and managers maintain the integrity, consistency, and quality, or multiple databases used in production, including Genomes on a Tree (GoaT), Sample Tracking System (STS), Collaborative Open Plant Omics (COPO), and BioSamples.

💻 **Bioinformatics**: Our bioinformaticians develop the suite of analysis pipelines that will run on every genome produced in Tree of Life, providing a central database of core results available for all.

🔩 **Systems**: We develop and maintain some core systems used in production, including the execution and tracking of all bioinformatics pipelines, and the deployment of third-party web applications for internal use.

### Tech stack {#it-tech-stack}

<img src="https://www.sanger.ac.uk/wp-content/uploads/Informatics-Infrastructure-Technology-wheel-300.png" alt="Tech stack" height=300px align="right">
The team uses a wide range of technologies, frameworks and programming languages, including Nextflow, Python, Conda, Jira, LSF, Singularity, and Kubernetes. The technology wheel below shows most of their logos. How many can you recognise?

### Members {#it-members}

Current members are listed below:

- [Matthieu Muffato](#matthieu-muffato)
- [Guoying Qi](#guoying-qi)
- [Cibin Sadasivan Baby](#cibin-sadasivan-baby)
- [Cibele Sotero-Caio](#cibele-sotero-caio)
- [Paul Davis](#paul-davis)
- [Beth Yates](#beth-yates)

#### Matthieu Muffato, Team Lead {#matthieu-muffato}

[<img class="emoji" src="/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/muffato-matthieu/) &nbsp; [<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/muffato) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/matthieu-muffato/) &nbsp;
</br> Matthieu leads the Informatics Infrastructure team, which guides the implementation and delivery of the genome assembly pipelines, and provides support for large-scale genome analyses for the Tree of Life faculty teams. He joined the Wellcome Sanger Institute in February 2021, to form the Informatics Infrastructure team for the Tree of Life programme. He has recruited 7 team members, with skills covering data curation & management, software development & operations, and bioinformatics.

#### Guoying Qi, DevOps Software Developer {#guoying-qi}

[<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/gq1) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/guoying-qi/) &nbsp;
</br> Guoying, a DevOps software engineer, has the responsibility of developing and deploying software and web applications for the Tree of Life project across various platforms such as computing farms, Kubernetes, OpenStack, and public clouds.

#### Cibin Sadasivan Baby, Senior Software Developer {#cibin-sadasivan-baby}

[<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/cibinsb) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/cibinsb/) &nbsp;
</br> Cibin, a Senior Software Developer, is tasked with designing and implementing the production systems for TOL-IT. Currently, Cibin is focused on building an automated platform to execute high-throughput genomic pipelines. The ultimate goal of this project is to develop a system capable of efficiently processing large amounts of genomic data.

#### Cibele Sotero-Caio, Genomic Data Curator {#cibele-sotero-caio}

[<img class="emoji" src="/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/sotero-caio-cibele) &nbsp; [<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/ccaio) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/cibele-sotero-caio-b379071a6/) &nbsp; [<img class="emoji" src="/assets/img/twitter.svg" height=25px>](https://twitter.com/CibeleCaio) &nbsp;
</br> Cibele is the data curator for the Genomes on a Tree (GoaT) - a platform developed to support the Tree of Life and other sequencing initiatives of the Earth Biogenome project (EBP).

#### Paul Davis, Data Manager {#paul-davis}

[<img class="emoji" src="/assets/img/logo/sanger-tol-logo-tree.png" height=25px>](https://www.sanger.ac.uk/person/davis-paul/) &nbsp; [<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/Paul-Davis) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/paul-davis-uk/) &nbsp; [<img class="emoji" src="/assets/img/twitter.svg" height=25px>](https://twitter.com/SirPaulDavis) &nbsp;
</br> Paul works on the main ToL Genome Engine. This system was developed by the ToL to manage and track samples from collection, onboarding, processing in the lab, sequencing and finally the publication of the assembly and Genome Note publication. As there are many steps in this process developing methodology to identify issues as early as possible is vital to avoid wasted time and resource. Paul works at all levels of the project fielding questions about data flow, data fixes and helps other ToL staff and project stakeholders with data and information. Paul also interacts with external groups and stakeholders to maintain data integrity in the public domain.

#### Beth Yates, Bioinformatics Engineer {#beth-yates}

[<img class="emoji" src="/assets/img/github.svg" height=25px>](https://github.com/BethYates) &nbsp; [<img class="emoji" src="/assets/img/linkedin.png" height=25px>](https://www.linkedin.com/in/bethanyates/) &nbsp;
</br> Beth is a Bioinformatics Engineer working on a building a platform to automate the production of [Genome Note publications](https://wellcomeopenresearch.org/treeoflife). The Universal Genome Note platform consists of a web portal, database and Nextflow pipelines. Beth is contributing to the genomenote pipeline, this pipeline fetches assembly meta data and generates some of the figures and statistics included in each genome note.

#### Alumni and Friends {#alumni}

- Alexander Ramos-Díaz, Google Summer of Code contributor. Alexander kick-started the development of the [BlobToolKit](/blobtoolkit) pipeline through the [2022 edition of Google Summer of Code](https://summerofcode.withgoogle.com/archive/2022/organizations/wellcome-sanger-institute).
- Zaynab Butt, Informatics and Digital Associate. Zaynab continued the development of the [BlobToolKit](/blobtoolkit) pipeline, and then developed the web interface of our internal tool for tracking workflow execution.
- Priyanka Surana, Senior Bioinformatician. Priyanka was a Senior Bioinformatician, overseeing the development of Nextflow pipelines for genome assembly, curation and downstream analyses. She is passionate about building networks that support peer learning and facilitated the workflows community on campus.

## Safety team {#safety}

The Safety team is responsible for ensuring the community is a safe place and responding to instances of misconduct. It is made up of community members who have displayed integrity, strong communication, and a genuine concern for community welfare.

The current safety team in 2023 are:

- [Priyanka Surana](mailto:ps22@sanger.ac.uk)
- [Matthieu Muffato](mailto:mm49@sanger.ac.uk)
- [Edward Symons](mailto:es13@sanger.ac.uk)

## Financial Support

<img src="/assets/img/contributors-colour/wellcome.png" alt="Wellcome" height=75px class="float-end darkmode-image me-5 mb-5 ms-3">
<img src="/assets/img/contributors-colour/sanger.svg" alt="Sanger" height=75px class="float-end darkmode-image me-5 mb-5 ms-3">

Tree of Life is a core funded programme within the Wellcome Sanger Institute and is allocated a slice of the Institute’s core Wellcome grant ref. 220540/Z/20/A.

<br clear="all" />

## History of sanger-tol

Established at the Wellcome Sanger Institute in 2019, the Tree of Life Programme is part of a global push to sequence the genomes of all eukaryotic life on Earth. That is, to produce reference-quality genomes for every animal, plant, fungi and single-celled protist living on our planet.

The initial focus of the programme was to provide whole genome sequencing capabilities for the Darwin Tree of Life project, a partnership of ten research institutions, with the aim of generating genomes for the estimated 70,000 species living in Britain and Ireland. Darwin Tree of Life makes a significant contribution to the Earth BioGenome Project, launched a year previously, an umbrella organisation of biodiversity genomics projects seeking to sequence all life. It also complements other like-minded projects, for example the Vertebrate Genomes Project which aims to sequence all vertebrates globally, and on which the Sanger Institute was already a partner prior to Tree of Life’s conception.

Almost half a decade and one pandemic later, Tree of Life has grown solid roots and many branches. The trunk of the programme is its genomics production pipeline, which supports collections of species in the field, their ethical and legal transportation to Sanger, extraction of high molecular weight DNA, and then the sequencing, assembly and curation of genomes organised at a chromosomal level to the highest quality. All this data is published freely and openly for researchers to use worldwide. With a demand to tackle so many different lifeforms, the innovation and expertise built by teams within the programme is vast.

Other projects have popped up alongside Darwin Tree of Life. For example, the Aquatic Symbiosis Genomics project which explores the genomics of symbiotic organisms in the ocean. Or BIOSCAN in the UK, which will study one million flying insects over five years to help better understand and monitor the health of their ecosystems. Tree of Life now has four full-time faculty, plus associate faculty members, exploring a multitude of genomic questions from rapid speciation to strange reproduction.

Key to being able to produce top-quality genomes at scale is building top-quality informatics systems. To support genome production, Tree of Life teams have either adopted existing softwares and tailored them to the programme’s needs, or developed bespoke platforms from scratch such as the Samples Tracking System. Within the bioinformatics teams, who piece together the DNA data into whole genome assemblies, new tools are being developed all the time to automate and increase efficiency, and to tackle specific issues, for example mitochondrial genome assemblies or contamination by DNA from other species. The breadth of talent in these teams ranges from early-career researchers to experienced hands whose time at Sanger can be traced back to the original Human Genome Project.

Together, these scientific endeavours, using the latest DNA sequencing technology, will enable us to tackle the extinction crisis, discover new biomedicines and biotechnologies, and better understand 3.5 billion years of life on Earth.
