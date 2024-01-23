---
title: 'Bytesize: nf-core/viralrecon'
subtitle: Sara Monzón and Sarai Varona - Instituto de Salud Carlos III, Madrid, Spain
type: talk
startDate: '2022-06-07'
startTime: '13:00+02:00'
endDate: '2022-06-07'
endTime: '13:30+02:00'
embedAt: 'viralrecon'
youtubeEmbed: https://www.youtube.com/watch?v=K1ThKn4p4u0
locationURL:
  - https://www.youtube.com/watch?v=K1ThKn4p4u0
  - https://doi.org/10.6084/m9.figshare.20020703.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/viralrecon

This week, Sara Monzón ([@saramonzon](https://github.com/saramonzon)) will talk about the newest developments in the nf-core/viralrecon pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=1)
(host) Hello, everyone. My name is Franziska Bonath. I'm the host of today's Bytesize Talk. And with me are Sarah Monzón and Sarai Varona from the Institute of Health, Carlos III. They are talking today about the nf-core pipeline viralrecon - updates and use cases. And up to you.

[0:28](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=28)
Thank you very much. Hello, everyone. We're really glad to be here today to talk about nf-core/viralrecon. This is the second talk regarding this pipeline, but the first one as a bytesize, I think. So we want to talk about some updates and some new functionality we've added to the pipeline in the last year and a half. And also, we want to speak about some use cases we've been using in our lab using viralrecon as the main character.

[1:04](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=64)
We're going to start with a brief road map, development map, followed by the major functionalities in the last releases. And as use cases, we're going to talk about the Relecov Network, which is the genomic surveillance network for SARS-CoV-2 in Spain, where we use viralrecon for data analysis. We're going to talk about a paper we participated regarding a study of a long-term COVID patient. And also, we're going to speak a little bit about the work we are currently doing, studying the multi-country monkeypox outbreak in ISCIII, where we are also using viralrecon.

[1:46](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=106)
As a roadmap, this pipeline was first started, the first release was in June 2020, but we started the development in March, more or less. The second major release was a year later, in May 2021, where the pipeline was rewritten using the DSL2 implementation. And also, a new whole branch of the pipeline was added for Nanopore data analysis. Also, Pangolin and nexclade was included for lineage assignment for SARS-CoV-2. Just a few months ago, in February, we created a new release, 2.3, that included some fixes regarding problems or decision-taking from ivar consensus. We're going to talk about this functionality deeper in the next slide.

[2:46](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=166)
Currently, we are in the version 2.4.1, and this is the major functionality we've added. The Nanopore branch of the pipeline allows us to handle both Illumina reads and Nanopore reads using viralrecon. For Nanopore data, the ARTIC Network pipeline is used. A variant calling and consensus genome output is generated and also Nextclade and Pangolin for lineage assignment is computed over this consensus genome.

[3:26](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=206)
One of the main functionality added in the version 2.3 is that now the user can determine which variant caller, which combination of variant caller or consensus generation software should be used. Until this version, by default, ivar variants always used ivar consensus as the software for generating the consensus, but now you can combine them. You can use ivar variants for variant calling and Bcftools consensus or the other way around, providing more flexibility for the user and also more capacity for the decision of how the consensus will be generated. This is one of the main functionality and this is important because it changes the output or the way the consensus is generated from previous versions.

[4:21](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=261)
Now the default is to use ivar variants as the variant caller and Bcftool consensus as the consensus generated. We've taken this decision due to some behavior of ivar that may not be desired for this case and some known issues of ivar consensus, that are not yet addressed by the software. For example, here we see that ivar includes low frequency selections and when we use viralrecon, we select a threshold for including variants in the consensus. For example, the default parameter is that we include variants in consensus that meet the criteria of more than 0.75 allele frequency. In this case, we see that even if we use this criteria in viralrecon, we see that this selection, which is 0.43 of allele frequency, is included in the consensus where the reference should be included. Here we can see the reference, the consensus generated by Bcftools and the consensus generated by ivar. Ivar is included and low frequency selected that shouldn't be there if you don't want to.

[5:42](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=342)
A known issue about ivar is that it has some issues regarding the calculation of default coverage of insertions and deletions. Here we can see that this is a low frequency deletion as the previous example. Again, the reference, the Bcftools consensus and the ivar consensus and here we see the low frequency deletion but an N and a mask position is added even if we have enough coverage in this area. So this is an issue about the depth of coverage calculation that is fixed using Bcftools consensus instead of ivar.

[6:27](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=387)
There is another issue why we selected Bcftools consensus. It is not an error of ivar consensus, it is just the behavior how ivar creates the consensus which may not be the behavior the user wants. That's the main difference between Bcftools consensus and ivar. If you want to include variants regarding intra-host via reliability or consensus, for example in this case, we have two positions here where ivar includes ambiguous nucleotides. This is because in this position, in order to meet the criteria of 0.75, ivar needs to add two nucleotides in this position. That's because we add the ambiguous nucleotide. In this case, if we only want to add the majority or the more representative nucleotide (in this case is A or G) these are the only two nucleotides that meet the criteria of more than 0.75. So it depends if you want to add all the information of intra-host via reliability in your consensus or you don't want to include this noise in your consensus. Ivar includes the ambiguous nucleotides because it includes the majority and the behavior is to include the majority alleles until you meet certain criteria. Bcftool consensus only includes variants that are more than a low frequency.

[8:15](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=495)
Another issue is just another example of the previous one. This is also an deletion in low frequency variants. We see that ivar is including ends, masking sequences that could make problems when you upload to GSAID. Instead of including the reference or the deletion. But this is an area well covered, but ivar is only including ends instead of the nucleotides or the deletion.

[8:56](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=536)
Another functionality we've added in this case, we are going to talk about two new functionalities regarding the script that converts the ivar output to VCF format. There we added two new functionalities regarding codon merging and strand bias. In the case of codon merging, we mean that when the variant of concern, B.1.1.7, appears in a new complex variant, it appears as a variant for SARS-CoV-2. And we realized that for this complex variant it changes the three nucleotides in a codon. The variant callers, ivar and all the variant callers reported the variants as three lines for three different changes. This is a problem because you don't have the correct annotation. These are three changes that change the codon entirely, so the amino acid is changed completely. If you have it in three lines, the annotation wouldn't be correct, not for ivar, not for SNPF. So we've created a function that goes position by position, reading the TCF file from ivar. And we check if they are consecutive. And if we found two or three positions that are consecutive, we check if they belong to the same codon. If they belong to the same codon, like this case, we see that the reference codon is exactly the same for the three positions. We collapse these three lines in just one. So the reference has the three alleles and the alternative has the three alleles. This makes that SNPF annotates the amino acid change correctly, fixing this issue. And this is included in the ivar variants script.

[10:56](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=656)
And another one is, as we all know, NGS data are prone for certain bias, strand bias is one of them. Here we found a strand bias. For example, when we have a variant that is only supported for forward or reverse strand reads. And this is normally creating a small probability that the variant is a false positive. A strand bias is usually corrected by most of the variant colors we have nowadays, but ivar still lacks this functionality. So we've added this annotation in the ivar output conversion to VCF. What we do is to create a contingency table regarding the forward and reverse strand reads for the reference and the right alleles. We calculate a first test and we mark as a strand bias position when the p-value is less than 0.05. This formula is obtained from the tutorial from GATK.

[12:07](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=727)
Finally, the new output for reporting variants is included also in the version 2.3. This is really useful because we combine the data from the variant calling, the anotation and also the lineage assignment. And this provides a good way to study, for example, metagenomics data from a sewage SARS-CoV-2 data. It is really useful for variant inspection, for studying co-infections, et cetera.
And now Sara is going to talk to you about the use cases.

[12:51](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=771)
Yes. Now I'm going to explain you three use cases of viralrecon here in the Institute of Health, Carlos III. The first one is the RELECOV Network, which is founded by the HERA-Incubator program and is a Spanish network that aims to create a SARS-CoV-2 surveillance at national level based on genomic sequencing.
In this network, the microbiology labs from hospitals are going to select the samples to be sequenced based on criteria established by public health authorities. And they are going to sequence those samples. Then they will send the FASTQ files to the RELECOV platform here in the ISC. And we are going to analyze those samples with viralrecon.

[13:43](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=823)
We are going to be able to see the national evolution of the viral variants and viral lineage. Also we are going to share genomic data with databases such as GISAID or ENA. The idea is that we will give support and information to the different labs that are inside the RELECOV Network.
As you can see in this schema, there is at least one group in each of the autonomous communities in Spain that are included in the network. So all together, we are going to create a national surveillance of SARS-CoV-2 and probably learn from this approach how to extend it to other pathogens.

[14:34](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=874)
This would be a general schema on how the samples are sequenced and analyzed here in the Institute of Health Carlos III, after two days of sequencing samples. They are going to be stored in a hard disk cabin and processed in a high-processing computing server here in the Institute of Health using viralrecon. And then the results are going to be redirected to the microbiology labs.

[15:05](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=905)
The second example is about an immuno-depressed woman that had prolonged viral replication. She was receiving immuno-chemotherapy. And after the last cycle of immuno-chemotherapy, six months after, she was admitted in the hospital after being positive for SARS-CoV-2 in a RT-PCR. After nine months of being discharged and re-admitted in the hospital, being RT-PCR positive for SARS-CoV-2 and receiving antiviral drugs and convalescent plasma, the woman died. What we saw after 257 days of collecting 12 samples for sequencing is that the last sample obtained had accumulated 29 nucleotide mutations and 22 amino acid mutations using viralrecon in the mapping approach with the Wuhan reference genome.

[16:11](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=971)
For this, we used viralrecon version 1.2 in development version. Something interesting is that we used the long variant table that Sara explained to create this plot where we selected the low frequency variants to see how they were changing over time in this patient. We have in the x-axis the date of sample collection and in the y-axis the allele frequency. Each line and dot represents one variant in the sample. When no dot is shown it means that that position didn't have enough coverage in the sample. So in this example, we can see that the ORF1AB mutations are mostly present in the non-structural protein 3. Something similar happens with the S-gene where most of the variants are accumulated in the spike protein S1. Also, we found one of the variants that was afterwards considered as a mutation of concern of the delta variant inside this woman when the delta variant wasn't circulating in Spain.

[17:38](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=1058)
Something interesting we found selecting this low frequency variants is that we saw patterns of different viral subpopulations competing intra-host. So we think that there was an intra-host mutation and competition between the virus subpopulation and also that those antiviral drugs were selecting resistant viruses.

[18:11](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=1091)
The last example is the most recent one and it describes how we at the Institute of Health Carlos III treated the multi-country monkeypox outbreak in non-endemic countries. We sequenced 28 samples and we used the latest version of viralrecon to obtain different FASTA genomes for both de-novo assembly and mapping approach against three different monkeypox genomes. We obtained, using an Illumina NovaSeq with 1x150 reads, 33 samples that had 100% of the reference genome covered, with at least a ten-fold depth. We used the mapping consensus FASTA files and de-novo assembly FASTA files to create multiple sequence alignments and see the performance of both approaches. We saw that the ends of the reference genome couldn't be assembled with the de-novo assembly approach, but with the mapping approach we could see that there was enough coverage to obtain those sequences.

[19:37](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=1177)
In the plasmid ID plots obtained with viralrecon we can see how in the de-novo assembly the right and left ends of the reference genome are missing. A monkeypox genome has short tandem repeats that we were trying to discover if this approach was able to obtain the exact number of short tandem repeats in our samples. We found that in the de-novo assembly approach when the short tandem repeats were inside different contigs, abacas introduces ends in between the context so we couldn't reconstruct the real short tandem repeat scaffold. In the mapping approach we saw that we had enough coverage to cover the reference short tandem repeats but that we are limited to the number of SDRs present in the reference. So in order to discover the real number of repeats present in our monkeypox samples we are trying to sequence the best covered sample with MiSeq 2x300. This is going to be analyzed with the latest version of viralrecon and with Oxford Nanopore technologies. We keep working on that so we can't tell you anything yet.

[21:06](https://www.youtube.com/watch?v=K1ThKn4p4u0&t=1266)
Well this is everything thank you very much for your attention. Thank you to all the people that developed viralrecon with us and to the reference laboratories in the institute and the economic unit for all this work. Thank you.
(host) Thank you very much. So now we have time for some questions for anyone. No? If there are no questions I also want to mention that you can always ask questions later on on the slack bytesize channel and this video will be uploaded to youtube later. Thank you very much again and I also would like to thank the Chan Zuckerberg Initiative for funding these talks.

</details>
