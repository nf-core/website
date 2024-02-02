---
title: 'Bytesize: variantcatalogue'
subtitle: Solenne Correard, University of British Columbia, Vancouver, Canada
type: talk
startDate: '2023-02-21'
startTime: '13:00+01:00'
endDate: '2023-02-21'
endTime: '13:30+01:00'
youtubeEmbed: https://www.youtube.com/watch?v=Em1cHCLQQ_c
locationURL:
  - https://doi.org/10.6084/m9.figshare.22140854.v1
  - https://www.youtube.com/watch?v=Em1cHCLQQ_c
---

This week, Solenne Correard ([@scorreard](https://github.com/scorreard)) is introducing us to the nextflow pipeline `variantcatalogue`. Though not an nf-core pipeline (yet), `variantcatalogue` is using nf-core derived material to aid in the creation of the pipeline.

The variant catalogue pipeline is a Nextflow workflow designed to generate variant catalogues (also called variant libraries): a list of variants and their frequencies in a population, from whole genome sequences. More information on the pipeline can be found [on Github](https://github.com/wassermanlab/Variant_catalogue_pipeline) or in the [associated preprint](https://www.biorxiv.org/content/10.1101/2022.10.03.508010v2).

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=1)
Hello, everyone, and welcome to this week's bytesize talk. The speaker today is Solenne Correard from the University of British Columbia in Canada, and she is going to talk about variantcatalog. This is a Nextflow pipeline, but it's not part of nf-core yet. Variantcatalog is used for population analysis from whole-genome sequencing, and specifically to identify variants and their frequencies. Since Solenne is living in Canada, and due to the big time difference, we decided that it's best to record this talk. Therefore, if you have any questions, please ask them in the Slack bytesize channel. As usual, I would like to thank Solenne for her time and the Chan Zuckerberg Initiative for funding the bytesize talk series. But now, without further ado, I hand over to Solenne.

[0:52](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=52)
Hi, everyone. Welcome to this week's nf-core bytesize talk. I am Solenne Correard. I'm a research associate at BC Genome Sciences Center in Vancouver, Canada. Today, I'm going to talk to you about the variantcatalogue pipeline. First, I would like to acknowledge the lands on which I work, live, and play. Those are the traditional, ancestral, and unceded territories of the Musqueam, Squamish, and Tsleil-Waututh nations.

[1:22](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=82)
First, what is a variant catalogue or a variant library? When we talk about genomics and DNA, a variant catalogue is the frequency of the variants within a population. For example, in this population of five individuals, they all get their whole genome sequenced and at a certain position in their DNA, some individuals carry an A, and some individuals carry a C. From that individual information, we can deduce the frequency of each allele in the population. In this example, the A allele has a frequency of 0.06, and the C allele has a frequency of 0.04. This is the main information that is within the variantcatalogue pipeline, the frequency of the variants within the population.

[2:08](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=128)
When do we use a variant catalog? There are several ways to use it, but a very good example is through the nf-core/raredisease pipeline. During the variant annotation and prioritization step of that pipeline, they use gnomAD. GnomAD is the biggest variant catalogue to date. The reason they use it is because a variant that is frequent in a population is unlikely to be responsible for a rare disease. When we are looking for the variant responsible for the rare disease in a kid, we can already filter out all the variants that are frequent in the population. As I mentioned, gnomAD is the biggest variant catalogue to date. It helped tons of families to get a diagnosis in rare diseases. But when we look at the ancestries of the individuals within gnomAD, we can see that most of the individuals are from European ancestry. Some populations are not even represented. They are not represented or underrepresented. This lack of representation from some population is leading to an inequity in genomic care. Because if the kids affected with a rare disease is from an ancestry that is not represented in the variant database, then it's harder to remove the variants that are frequent from this population, and so harder to give a diagnosis to this kid.

[3:47](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=227)
This is a known issue. Several variant catalogs were generated around the world. For example, Iranome with the Iranian population, or KOVA with the Korean population. The project I was working on is the Silent Genomes Project in Canada. It's a partnership with the Indigenous populations of Canada to build the Indigenous background variant library. A very similar project is taking place in New Zealand with genomics Aotearoa, where they are working with the Māori population. When we were working on the Indigenous background variant library, we needed a pipeline to process the data to get the variant frequencies.

[4:30](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=270)
Some pipelines existed, but none of them were fulfilling the three constraints that we had. The first one is that we wanted the pipeline to rely on open access tools that were previously benchmarked, because we didn't want to develop any new software or new tool. We wanted it to be comprehensive. By that, I mean it had to include single nucleotide variants, but also mitochondrial variants, structural variants, short tandem repeats, and mobile element insertions. All those classes of variants are known to be potentially implicated in rare disease. It's very important that all of them are present in a variant catalogue. Finally, we wanted it to be able to work on local servers or on the cloud, because different projects may have different constraints.

[5:23](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=323)
We developed the variantcatalogue pipeline that you can see on the left here. This is just an overview, and I'm going to describe each part in more details. But the idea is that it takes as input FASTQ files from participants. It outputs VCF files, so variant calling files with information about the variants, the position, the allele, the frequency of this variant within the population, which is the main information we want, the frequency by sex, as well as some annotation. The pipeline is divided in four subworkflows that can work independently, or all of them can be run in parallel, or at least in the same pipeline.

[6:09](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=369)
The first subworkflow is the mapping subworkflow. It takes as input short read paired-end sequences for individuals, as well as a reference genome. It has been developed so far for GRCh37 and GRCh38. The mapping tool is bwa mem, and it outputs one BAM file per individual. The second subworkflow is the mitochondrial variant subworkflow. It is very much based on the work from Laricchia et al. that was published in 2022. It's therefore very similar to the pipeline that is used by gnomAD for their mitochondrial variants. It takes as input the BAM files previously generated, the variant caller for the mitochondrial variants is GATK Mutect2.

[7:06](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=426)
The reason why there is a parallel section here is because the mitochondrial DNA is circular. To be able to map the mitochondrial reads against this reference genome, it is linearized with a fake breakpoint around zero here. The reads that are supposed to map over the fake breakpoint do not map correctly, hence the variants located around these regions are not called correctly. To address that issue, they developed a shifted reference genome where the fake breakpoint is located on the other side, which allowed to call the variants correctly in this region. Then the variants are lifted over, the information is merged into several VCFs. I will detail the steps at the bottom later.

[8:04](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=484)
The third subworkflow is the single nucleotide variant subworkflow. It is the most straightforward one. For variant calling, we decided to use DeepVariant. We are using GLnexus for the joint calling. For the fourth subworkflow, which is the structural variant subworkflow, it was mostly developed by Mohammed Abdallah, a postdoc within the Wasserman Lab. It was decided to use Smoove and Manta for structural variant callers. Jasmine is used to merge the variants, and Paragraph is used to genotype the structural variants within the individual data. Then the information is merged with BCFtools. For the short tandem repeats, we are using ExpensionHunter. For the mobile element insertions, we are using MELT.

[8:57](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=537)
All the variant calling part is very similar to other pipelines, such as nf-core/raredisease or nf-core/sarek. What is really specific about this pipeline is the steps at the bottom here. It's the sample quality control, variant quality control, allele frequency calculation, and also sex imputation. The reason for that is the quality control is performed differently if you have just one individual or a trio versus if you have a population. All of this is performed within Hail, which is a Python-based analysis tool that is also used by gnomAD and some other variant catalogue pipelines. As I said, it performed some quality control as well as the variant frequency calculation. Then the variants are annotated using VEP.

[10:00](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=600)
That was just an overview of the pipeline. This is the actual complete pipeline. It's available on the Wasserman lab GitHub, and it's described in more detail in this preprint. It was tested on 100 samples, and it works. The details of the number of CPU hours, as well as the number of samples and variants that were filtered out by the quality control steps is available within the preprint. However, this version still rely on locally installed software. That is an issue for two reasons. First, it's really hard for other projects to use. Second, it's impossible to test very easily. We are used to test other nf-core and Nextflow pipelines with just one command line.

[10:48](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=648)
Therefore, the future for the pipeline is to move it to an nf-core level pipeline. My goal is to move the mapping as well as the single nucleotide variants workflows during next month's hackathon. If anyone wants to team up with me for the code or for coding review, please reach out. After that, we will have to move the mitochondrial and the structural variants workflows also to nf-core. This will allow first other people to try it more easily, but it will also force us to do better documentation. That is very important to make sure that other groups can use the pipeline. If the documentation is good, then it's easier for other people to try and use this pipeline.

[11:40](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=700)
To test the pipeline, I actually needed to create a new data set because the one that were available within nf-core did not fit my needs. I needed paired-end short read FASTQ files that included part of an autosome as well as parts of chromosome X and Y to impute sex for the individuals as well as read mapping to the mitochondrial chromosomes to test subworkflow 2. I also needed reads supporting the presence of a structural variance to be able to test the subworkflow 4. And several samples, including XX and XY individuals, to be able to test the variant frequency calculation part. This will hopefully be available to others soon, in case you need them to test your tools or your pipeline. I will also include the reference genome for the same region and additional files such as the short tandem repeat catalog, the mitochondrial reference file and the shifted one I mentioned before.

[12:46](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=766)
In other future developments, I would like to include more reference genomes, including the T2T for humans but also non-humans reference genomes. I would like to include more software, for example, to give the opportunity for the user to decide which mapper they want to use, which variant callers they want to use. We also want to make sure that it fits with the nf-core/raredisease pipeline. I know that they use slightly different callers for structural variants. It would be interesting to make sure that there is a good fit. It's also possible to include additional metrics such as ancestry inference, mitochondrial group assignment or relatedness calculation. Those are metrics that are often associated to variant catalogue pipelines. It was out of the scope for the Silent Genomes Project, but we understand the relevance for other projects and it would be great to also include them and have them as an option.

[13:53](https://www.youtube.com/watch?v=Em1cHCLQQ_c&t=833)
I would like to acknowledge everyone within the Wasserman Lab, especially Wyeth Wasserman, the team leader, Mohammed Abdallah, who worked a lot on the structural variant pipeline subworkflow and the rest of the pipeline, as well as Brittany Hewitson, the Silent Genomes team, and also all the nf-core community. It's been a very welcoming community and I've learned a lot. Obviously, this is not live. If you have any questions, please reach out on the nf-core/variantcatalogue channel to sparkle a discussion and start threads on different things. If you prefer to reach out directly to me, you can do it through Twitter or GitHub. Thank you for your attention and have a great rest of your day.

</details>
