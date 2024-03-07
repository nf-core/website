---
title: 'Nextflow in Action'
subtitle: Zahra Waheed and Ahmad Zyoud; Jonathan Manning
type: talk
start_date: '2023-09-05'
start_time: '15:00+01:00'
end_date: '2023-09-05'
end_time: '16:00+01:00'
location_name: BIC Conference Room
---

> Organised by the [the BioDev Network](https://www.youtube.com/@biodev-network)

# Nextflow in Action

Join us for an engaging exploration of Nextflow in designing, executing, and managing complex bioinformatics workflows. Our meetups are a platform to exchange insights and delve into practical use cases. If you missed the previous sessions, you can watch the recordings on [YouTube](https://www.youtube.com/playlist?list=PLo5QmrytFHLHUkBLviJykEHYE8ZKzJOm5). This session includes:

##  Zahra Waheed and Ahmad Zyoud – Data Coordination and Archiving, EMBL–EBI

The COVID-19 pandemic highlighted the importance of sharing genomic data and metadata globally, through submission to public nucleotide sequence repositories. To date, over 12 million raw reads and assembled SARS-CoV-2 sequences have been submitted to the European Nucleotide Archive (ENA) alone, which are visible and retrievable from the COVID-19 Data Portal. The rapid sharing of this data (whilst keeping in line with recommended metadata standards) is key to efficient outbreak surveillance and sequence interpretation and helps to drive a more effective public health response. Thus, developing simple, user-friendly submission tools is valuable for lowering the barrier to data entry, and maximising the rate and volume of data shared for scientific research.

Seeing the need for a one-stop shop SARS-CoV-2 submission tool, the ENA team developed the SARS-CoV-2 Drag and Drop Uploader, which requires no technical skills from users and no prior knowledge of the repository’s submission process. The tool offers an alternative route to submit SARS-CoV-2 data to the ENA but can also be repurposed for other viral submissions, as proven by the equivalent ENA Monkeypox Uploader in response to the 2022 outbreak.

Here we present a full Nextflow pipeline for the back-end automation of the SARS-CoV-2 Drop Uploader, incorporating existing ENA APIs, command line programs and AWS. We additionally present a portable, standalone version of this workflow for general ENA data submission.

<iframe width="560" height="315" src="https://www.youtube.com/embed/gm06vgrgijc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

## Jonathan Manning – nf-core

[nf-core/differentialabundance](https://nf-co.re/differentialabundance) is a bioinformatics pipeline that can be used to analyse data represented as matrices, comparing groups of observations to generate differential statistics and downstream analyses. The initial feature set is built around RNA-seq, but we anticipate rapid expansion to include other platforms.

(No recording available)
