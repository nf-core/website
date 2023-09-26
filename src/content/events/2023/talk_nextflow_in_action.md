---
title: 'Nextflow in Action'
subtitle: Tree of Life, Wellcome Sanger Institute and Genomics England
type: talk
start_date: '2023-09-26'
start_time: '14:00+01:00'
end_date: '2023-09-26'
end_time: '15:00+01:00'
location_url: https://sanger.zoom.us/j/96537565366?pwd=UWpYVVdsbTZFREZkbndaZGNkNGJrdz09
---

# Nextflow in Action

Join us for our **monthly series** of talks: **“Nextflow in Action”**.

**20 minutes** for each talk + questions, we will be focussing on topics about how different organisations are using and developing with Nextflow.
These will be recorded and made available on [YouTube](https://www.youtube.com/@workflows-community).
It is our hope that these talks / videos will help us learn from the greater Nextflow community. If you want to present, please reach out to [Priyanka Surana](mailto:ps22@sanger.ac.uk).

## Matthieu Muffato, Tree of Life, Wellcome Sanger Institute
**Adopting Nextflow in Sanger Tree of Life**
<br>Tree of Life is the latest scientific programme of the Wellcome Sanger Institute, aiming at producing reference genome assemblies for tens of thousands of species. We decided early on that Nextflow would be the workflow manager of choice - what an adventure it has been ! In this presentation, I will tell you about our journey, starting from training users, setting up development environments, right up to advertising our complete pipelines and reaching out to users outside of the department, with some insights into deploying and automating the pipelines.

## Luke Paul Buttigieg and Ricardo Humberto Ramirez Gonzalez, Genomics England
**Porting clinical genome analysis to Nextflow at Genomics England**
<br>Genomics England provides whole genome sequencing diagnostics to the Genomic Medicine Service (U.K), a free at the point-of-care, nationwide, genomic diagnostic testing service, with ambitious targets of processing 300,000 samples by 2025. Currently, all clinical bioinformatics is processed using a clinical-standard certified, internally developed workflow engine (Bertha). We are migrating to a new solution (Genie) which combines Nextflow and Nextflow Tower with custom functionality, so we can focus on our core mission to enable equitably accessed, genomics medicine for all. Genie should help us support newer use cases quicker, across different infrastructures such as cloud, and uses a standard workflow definition language. We have developed an approach to migrate at speed in an agile and iterative fashion. We are mocking the services in the workflow environment in Docker containers to allow us to run continuous integration tests. We are using an automated comparison testing framework to compare the existing system with the new one to detect regressions. Later, we will iteratively refactor the workflows, breaking up the Singularity image and optimising for performance. In this talk, we will describe this migration strategy, risk management and lessons learnt while working through this large-scale effort.

