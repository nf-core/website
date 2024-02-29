---
title: 'Nextflow in Action'
subtitle: Matthieu Muffato, Luke Paul Buttigieg and Ricardo Ramirez
type: talk
startDate: '2023-10-12'
startTime: '14:00+01:00'
endDate: '2023-10-12'
endTime: '15:00+01:00'
location_name: Green Room, Hinxton Hall, Wellcome Genome Campus
location: https://sanger.zoom.us/j/96537565366?pwd=UWpYVVdsbTZFREZkbndaZGNkNGJrdz09
youtubeEmbed: https://youtu.be/w4RGam3T8iQ
---

**Matthieu Muffato, Tree of Life, Wellcome Sanger Institute**

_Adopting Nextflow in Sanger Tree of Life_

Tree of Life is the latest scientific programme of the Wellcome Sanger Institute, aiming at producing reference genome assemblies for tens of thousands of species. We decided early on that Nextflow would be the workflow manager of choice - what an adventure it has been ! In this presentation, I will tell you about our journey, starting from training users, setting up development environments, right up to advertising our complete pipelines and reaching out to users outside of the department, with some insights into deploying and automating the pipelines.

**Luke Paul Buttigieg and Ricardo Humberto Ramirez Gonzalez, Genomics England**

_Porting clinical genome analysis to Nextflow at Genomics England_

Genomics England provides whole genome sequencing diagnostics to the Genomic Medicine Service (U.K), a free at the point-of-care, nationwide, genomic diagnostic testing service, with ambitious targets of processing 300,000 samples by 2025. Currently, all clinical bioinformatics is processed using a clinical-standard certified, internally developed workflow engine (Bertha). We are migrating to a new solution (Genie) which combines Nextflow and Nextflow Tower with custom functionality, so we can focus on our core mission to enable equitably accessed, genomics medicine for all. Genie should help us support newer use cases quicker, across different infrastructures such as cloud, and uses a standard workflow definition language. We have developed an approach to migrate at speed in an agile and iterative fashion. We are mocking the services in the workflow environment in Docker containers to allow us to run continuous integration tests. We are using an automated comparison testing framework to compare the existing system with the new one to detect regressions. Later, we will iteratively refactor the workflows, breaking up the Singularity image and optimising for performance. In this talk, we will describe this migration strategy, risk management and lessons learnt while working through this large-scale effort.

