---
title: 'Bytesize: nf-core/funcscan'
subtitle: Jasmin Frangenberg - Dept. Palaeobiotechnology, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute
type: talk
startDate: '2023-01-24'
startTime: '13:00+01:00'
endDate: '2023-01-24'
endTime: '13:30+01:00'
embedAt: 'funcscan'
youtubeEmbed: https://www.youtube.com/watch?v=c1CnE6jPhpg
locationURL:
  - https://www.youtube.com/watch?v=c1CnE6jPhpg
  - https://doi.org/10.6084/m9.figshare.21953978.v1
---

This week, Jasmin Frangenberg ([@jasmezz](https://github.com/jasmezz)) is going to introduce nf-core/funcscan. nf-core/funcscan is a bioinformatics best-practice analysis pipeline for the screening of functional components of nucleotide sequences such as assembled contigs. This includes mining for antimicrobial peptides, antibiotic resistance genes and biosynthetic gene clusters.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=1)
(host) Hello everyone and welcome to this week's bytesize talk. With us is Jasmin Frangenberg. I'm very happy that you're here. Thank you very much. She's going to talk about yet another new pipeline that is going to be released very soon, which is nf-core/funcscan. Off to you, Jasmin.

[0:22](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=22)
Yes, thank you very much. I will introduce this pipeline to you now, which is an nf-core pipeline to screen for functional components of nucleotide sequences from prokaryotic genomes or metagenomes. What are these functional components that we are interested in or that we screen for? The pipeline screens on the one hand for antimicrobial peptides. These are important in innate immunity and they are very short sequences, peptides out of about 20 amino acids, so you can find them even in small or fragmented DNA and metagenomes. The same applies to antibiotic resistance genes. On the other hand, biosynthetic gene clusters, here at the bottom. They are quite big, because they consist of a whole gene cassette, which codes for a whole metabolic function, secondary metabolites or natural products.

[1:24](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=84)
Who would be interested in such a pipeline, which identifies these compounds? In natural product discovery, where you can identify these compounds to develop therapeutics, in antibiotic research, in environmental metagenomics, or simply to have functional and genomic annotations. In these research fields, the detection of these compounds is already being done with a couple of tools, however, there are certain issues. One of them would be the efficiency, because mostly you apply the tools manually and then you only have a very specific purpose of the tool. You can identify a single compound, but it's not very broad and you have only a single algorithm that identifies the output. It could be more feasible to have this whole process streamlined in a pipeline. Also the output of these tools is not standardized. Another issue would be the reproducibility, because throughout the years, the tools develop new functions, bugs are fixed. It's very important for researchers to record which versions of which tools they are using, which is hard if you execute them manually on your samples. Also data privacy, there are a bunch of tools that offer web services where you can upload your data where they are analyzed for you. However, this requires that you give your data to a third party, which is not always intended or even possible. Another issue is that bioinformatics skills are often needed. Sometimes you even have to write small bash scripts to execute the tools on your data, which is not possible for all people. For example if they are biochemists who just want to know what is in the data, they don't want to be trained bioinformaticians.

[3:30](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=210)
These are the many problems that our pipeline tackles, namely that it is very scalable since it's a Nextflow pipeline, all nf-core pipelines are Nextflow pipelines. They are very efficient and scalable. You can execute them on your local computer, laptop, up to the institute's HPC. They are reproducible since they record all the tools and versions of the tools. Of course, you can decide where you want to have your data, you are not forced to put them on any web server. Also, it is very easy to execute the pipeline, which you will see later when we come to the tutorial part.

[4:11](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=251)
I emphasized how easy the pipeline is to use, but it didn't start very easily. I go back to October 2021, when we assembled the ideas to develop a pipeline of many tools, we brainstormed what would be needed for obtaining the resistance genes, the biosynthetic gene clusters and the AMPs. Not all tools were yet on Conda or had nf-core modules. We had to do a lot of work there. Then throughout the next year, we streamlined the process a bit and the ideas got clearer. We even made the first sketch of the famous tube map sketch. Finally in 2023, the pipeline is ready to use. This is the current version.

[5:05](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=305)
I will walk you through it. In the first step, we have the input which is being annotated. As I said, input can be any genome sequence, could be metagenome, contigs, could also be complete bacterial genomes. This data is then analyzed by one of the three tools, the annotation tools. After this, this data goes into one or all of the three workflows. The antibiotic resistance genes in the yellow workflow, the BGCs in purple and the antimicrobial peptides in red. Not all of the downstream tools need the annotated data. For some, we also use the direct input data.

[5:50](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=350)
Then as I said, each of the workflows has a bunch of tools. For example, the AMP workflow has four tools. As I mentioned before, they follow different strategies. Some of them use, for example, deep neural networks and machine learning to identify compounds of AMPs, which would be, for example, ampir or here DeepBGC for the BGC workflow. Other tools have rule-based strategies. There are a lot of algorithms predicting the compounds and the results are then very diverse, as you can imagine. It is important to aggregate these outputs and summarize them into a nicely readable format which is the third step. For this, we use one tool per workflow, two of them are developed by ourselves – AMPcombi and comBGC – and hAMRonization was already a tool available.

[6:51](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=411)
This was the overall workflow and now I would like to show you how to apply the pipeline and you will see that is really very easy. We start with the input, which is a sample sheet, basically a table with two columns. The first one is your sample name, the second one is the path to your FASTA input file. Of course your FASTA file includes the ID of your sequence and the sequence itself. This is what you need to actually run the pipeline and it is as easy as running `nextflow run nf-core/funcscan`. You give your input sample sheet, give your output directory. This is a minimal example of a pipeline run. Of course, it is recommended to use more parameters. One of them would be in the annotation step, the flag `--annotation_tool`, where you can decide which tool you want to use. They have different properties. For example, prodigal is very fast, however, we noticed that with prokka we get better downstream results. It depends on your needs and ideas, which tool you would like to choose. The default is prokka.

[8:00](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=480)
After the annotation step, we come to the actual identification of the compounds. You can activate each workflow with this flag `--run_amp_screening`, for example, for the AMPs. And by activating this, all the AMP tools are run on your data. You can also choose, for any reason, to deactivate any of the tools. You can switch them off with the flag `--amp_skip` and then the name of the tool. This might be because some tools might be very slow or you think they are so specific that you are not interested in the output. As I said, for whichever reason, you can switch them off. This is the same for the antibiotic resistance workflow. You can apply this flag, it runs all the four or five tools on your data and you can skip any tool with the `--arg_skip` flag. Same applies for BGC identification. You have the flag, all the tools are run, you can skip whichever you might want to skip. Of course, you can use not only one of the flags per run, but all three flags at the same time. Your data is investigated simultaneously and parallelized as much as possible with Nextflow. Okay, so these are the identification steps.

[9:21](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=561)
Now we come to the summary steps for each workflow. Let's start with the antibiotic resistance, which is done by hAMRonization, which is a tool that is already out there. Here you can see the GitHub link. This tool can actually summarize a bunch of outputs of resistance identification tools. Our pipeline currently includes the orange tagged ones. The output of those tools is then summarized into a standardized gene report. This is how it looks. It's a table with a lot of columns. You have here the sample IDs, then the genes that have been identified, some information about the databases, which tools were run, and so on. These are actually all the column headers that are very conclusive and you can use this output table for downstream analysis in R or any statistics program.

[10:17](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=617)
This is very similar to AMPcombi, which we developed ourselves, Anan and Louisa developed this, where you also have your sample IDs and then some information about probability of AMPs. Additional feature is that it not only identifies your antimicrobial peptides, but it also does some back aligning to a reference database to identify taxonomic classification. It also infers some chemical properties like stereochemistry and provides the publication so you can go back and read more about the compound identified. The last tool for the BGC workflow is comBGC. Similar fashion, we have the sample IDs, the tools which have been applied, and then more information about your candidate biosynthetic gene clusters. With this, you see that we have a scalable workflow now to identify these compounds, which are important for a couple of research fields for, as I said, drug development, antibiotic research and so on.

[11:28](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=688)
Since the pipeline is almost ready, it's probably going to be released next week. Let's see about it. We have at least added all the modules and subworkflows. We do some more testing and then the pull request will go out. I can already advertise if there is someone here in the chat, who would like to review, please feel free to reach out to us on Slack. In the future, we would like to include more screening modules and to also have a visual summary of the output, which would be a graphical dashboard, probably with a Shiny app. Let's see about that.

[12:11](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=731)
With that, I would like to introduce the development team, which is James, Louisa, Anan, Moritz and me. Of course, we got a lot of help from the nf-core community, which was always assisting, very nice community. Also I would like to emphasize some colleagues here at my institute, which helped with biological and biochemistry knowledge. My supervisor, Pierre Stallforth from the Leibniz HKI. With this, I would like to close and lead you to our repository and the documentation of the pipeline. If you want to interact with us, feel free to join us on Slack and otherwise I'm open for questions either now or later on Slack. Back to you, Franziska.

[13:03](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=783)
(host) Thank you very much. Very interesting. Anyone can now unmute themselves if they have any questions, they can also post questions in the chat and then I will read them out. Are there any questions from the audience? Otherwise I actually have a question.

(question) You have shown a minimal command that you can run, that doesn't actually specify the workflow that it's using. Is that going to use all three workflows or a specific one, a default?

(answer) This one you mean? Exactly. In the default we have specified none. This would actually run only the annotation, which is probably not very useful for you. This is the current state of the settings. Maybe we will change this later. I don't know.

[14:05](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=845)
(question) Right. Would it make sense to run all three workflows at the same time or is that different kinds of samples?

(answer) No, no, that's what it's designed for, to run efficiently on all three workflows. It depends on your interest: If you are not interested in the resistance genes, then of course you don't need to run it, but it's very efficient to use this also.

[14:26](https://www.youtube.com/watch?v=c1CnE6jPhpg&t=866)
(host) Thank you. Are there any more questions at this moment in time? Otherwise, I thank you again. It was a very nice talk. Of course I would also like to thank the Chan Zuckerberg Initiative for funding our bytesize talks and our audience for listening to the talk. I hope to see everyone next week. Thank you very much. Bye.

</details>
