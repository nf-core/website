---
title: 'Bytesize: nf-core/quantms'
subtitle: Julianus Pfeuffer (Zuse Institute Berlin) and Yasset Perez Riverol (EMBL-EBI)
type: talk
startDate: '2023-05-30'
startTime: '13:00+02:00'
endDate: '2023-05-30'
endTime: '13:30+02:00'
youtubeEmbed: https://www.youtube.com/watch?v=pBzelkgrPgQ
embedAt: quantms
locationURL:
  - https://www.youtube.com/watch?v=pBzelkgrPgQ
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/quantms

This week, Yasset and Julianus will introduce to us the nf-core pipeline quantms. nf-core/quantms is a bioinformatics best-practice analysis pipeline for Quantitative Mass Spectrometry (MS). Currently, the workflow supports three major MS-based analytical methods: Data-dependant acquisition (DDA), either (i) label-free or (ii) isobarically labelled quantitation (e.g. TMT, iTRAQ), and (iii) Data-independent acquisition (DIA) with label-free quantification.
More information about nf-core/quantms can be found on the [nf-core website](https://nf-co.re/quantms)

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1)
(host) Hello everyone, sorry for being late, but I'm very happy today that we have Julianus and sorry, I cannot pronounce your name.

(speaker 1) That's fine.

(speaker 2) Yasset.

(host) Yasset. Yes. To have a bytesize talk today, I will hand it over to you to introduce yourself if you want to.

[0:32](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=32)
(speaker 2) So I am Julianus, Julianus Pfeuffer. I am a postdoc at the Zuse Institute Berlin, and I was working a long time doing my Ph.D at the University of Tübingen in the Freie Universität Berlin. I was a long time OpenMS contributor and maintainer. A lot of the pipeline that we will present today is, of course, based on OpenMS and it will be about mass spectrometry.

(speaker 1) Okay, I can present also myself. I am Yasset Perez Riverol, I am the team coordinator of the Freie Database, which is the largest proteomics database of NBLEDM and I've worked with also one of the developers of quantms. If there are no further questions in the beginning, I also can start sharing my screen. Would that be okay? Franziska?

(host) Your co-host. You should be.

(speaker 2) Okay. Let's see. Do you see the presentation or do you see some of my PowerPoint things?

(host) We see the presentation.

(speaker 2) Perfect.

[2:30](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=150)
The presentation today will be about a Nextflow workflow this time and it's not a super recent addition to the nf-core community platform, but it was there for a long time or a little bit longer time, but we just recently released a 1.1 version that is much more stable and much more up to nf-core standards. We thought that would be a great time to introduce it. It is a workflow for, as the name implies, quantitative mass spectrometry data analysis. Of course, following all the nf-core standards, it is meant to be very reproducible and also applicable for large scale analysis. For example, on those big public repositories like PRIDE, where Yasset is from.

[3:43](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=223)
What are the areas of application of our workflow? It is planned as a most-in-one workflow for the analysis of quantitative mass spectrometry experiments in general, that can mean metabolomics, proteomics, proteogenomics, but the current focus and the topic that we started with in the workflow is relative quantification of proteins or modified proteoforms based on mass spectrometry experiments. I will present first or ask the question, why do people do proteomics and how is it different from the usual genomics that we see here in nf-core? One nice example that people always give is the difference between a caterpillar and a butterfly. While they share exactly the same genome, unless some very slight mutations due to the stage of life, they have a vastly different proteome, not only in the amount of proteins that are expressed, but also the types of proteins that are expressed and how they are modified. All of this gives a much better representation of the actual phenotype of an organism or a cell. The proteomics can be used in addition or instead of genomics or even transcriptomics.

[5:26](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=326)
One technique to get the quantities of proteins, for example, in the cell is mass spectrometry, and I would say the most common technique to do that is via liquid chromatography coupled mass spectrometry. This means you digest your sample with the proteins first, you put it into an Eppendorf tube. You subject this to a liquid chromatography to split them up by some physical chemical properties, to make the analysis easier. You ionize them, put them into the mass spectrometer, and the mass spectrometer can then measure the ions and the amount of ions that are there based on their property or their behavior in a magnetic field. What you get out is a so-called mass spectrum, where you can see the intensity which is related to the amount of ions that were there for a specific mass or to be specific, a specific mass-to-charge ratio. The problem with this is that those mass spectrometry experimental setups can be very complex. People add new types of mass spectrometers, new types of experiments that they want to. They invent new labeling strategies and how to compare different samples and so on.

[7:16](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=436)
Here I give a little overview over the most common strategies in quantitative proteomics and want to highlight which of them are supported by quantms. Overall, we only support relative quantification since absolute quantification can usually only be done with certain standards and few people do it, unless they have a really exciting project going on. But also in relative quantification, you have two big subtypes. You can have labeled relative quantification and label-free relative quantification. Label-free is usually cheaper because the labels are expensive, but the analysis is sometimes a bit more complex. Here you can have so-called data-independent acquisition, data-dependent acquisition. You can have it feature-based or spectral counting. Regarding labeling, you can label the proteins or the protein pieces, the peptides, in vitro or you feed your organism or cells certain labeled amino acids in vivo. Our focus for labeled quantification was the so-called TMT and iDRAC strategies, which are very similar in analysis types.

[8:53](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=533)
If you have any data set that was gathered with one of those green strategies here, quantms should be very useful for you. Here's an overview of the pipeline. Everything starts with the spectra in mzML or in raw format, mzML is open. It's sometimes a bit more verbose, bigger. It's an XML-based format. But we can also read raw files from the Thermofisher instruments through conversion. We do some pre-processing on the spectra as well as on the protein database that you give the pipeline to say which proteins you think that are in your sample and you would like to identify and quantify them. Then we have three different branches, that depend on which strategy the experiment was based on. We have this data-dependent label-free branch in blue, data-dependent isobaric labeling, TMT or iDRAC in red, and data-independent acquisition in green. The top ones are usually done by OpenMS tools. It's a framework for mass spec analysis, while the lower one is done by DIN. It's a separate package, where we were in close collaboration with the author to make it as efficient as possible in a distributed computing environment.

[10:48](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=648)
I will go over the steps one by one. At first, a bit more details about the input. I mentioned the mass spectra already. The second one is an experimental design that you need, and we highly recommend to use the sample-to-data relationship format. It's a community-developed and tab-separated format for the data sets, for example, that are in PRIDE, and we annotate a lot of them manually for our re-analysis. It contains information about the contents of the samples, like organism, labeled or not, the experimental setup, but also the biological question, like which condition a sample belongs to. The protein database is the usual easy FASTA format, and can be either directly downloaded from SwissProd or trEMBL, or manually created by some proteogenomics studies you did before. It can be with or without so-called decoy proteins that we need later for false discovery rate estimation.

[12:08](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=728)
In the pre-processing, we, as I said, convert and index all of our spectra, and the default format will be our mzML. Everything else you have will be converted into mzML before, keep that in mind. We combine information from the SDRF and the Nextflow parameters. Currently it's a bit mixed where you can set certain parameters, because we also wanted to support very simple designs, where a lot of the information is implied. We do some sanity checks, convert them into designs for the specific tools, but also units for specific tools or certain vocabulary for specific tools. Regarding the database, we can also generate the decoys for you, it's usually done by reversing or shuffling sequences in the database. Then we perform identification with common so-called database search engines. Currently you can select between MSGF+ and Comet, or both of them, in which case they will be combined probabilistically by an OpenMS tool called ConsensusID. We then offer a rescoring mechanism that uses more features than just the similarity between predicted spectrum and observed spectrum. This is currently only possible by the SVM-based tool Percolator, but we are heavily developing or trying to integrate deep learning-based scores from, for example, MS2rescore. The false discovery rate estimation is done based on the well-established target decoy approach. We offer FDRs on multiple levels, the peptide-spectrum-match level, peptide-level, protein level, protein-group level, and on different scales, either for a specific sample only or for the whole experiment. We can do the so-called picked FDRs that were recently published and show a bit more sensitivity in large-scale experiments.

[14:42](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=882)
For the quantification of the peptides, in label-free quantification we use the OpenMS ProteomicsLFQ tool, which is also the main part of the old nf-core ProteomicsLFQ pipeline, which means if you're using that one, this is fully integrated and superseded, so you may switch to quantms. This performs the following tasks. It does the identification of quantifiable features in your mass spec data. This can be done targeted by looking for specific IDs or untargeted by just looking at isotopes and illusion shapes. It then does retention time alignment. It links the identifications to get the best matches over all samples, and then you can optionally also transfer identifications to features that do not have an identification, or you can re-quantify parts of your MS experiment if in all samples but this one, or in the most samples but this sample, there was a feature but you couldn't find one in this one, then you can extract the last part of the signal.

[16:13](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=973)
Then isobaric labeling, it's much easier because it's just based on the intensity of so-called reporter ions. We support most TMT and iTraq plexes, which means the plex just tells you how many channels you can multiplex into one sample, which means how many samples you can have in one mass spec run, let's say. We also support so-called SPS, which introduces a third fragmentation level for mass specs.

[16:56](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1016)
When you have quantified the peptides, you usually are interested in the proteins that they come from, therefore we have two different inference techniques implemented, the Bayesian one with the OpenMS tool EPIFANY, but also a simple rule-based aggregation of peptides to proteins. Regarding quantification, we support the common strategies, TOP3 peptides per protein. ibaq is a common strategy that normalizes by the length of the protein, for example. Those come from OpenMS, but we also have support for statistical post-processing tools like MSstats and Triqler, which then they have much more elaborate statistical models, and they also include significant testing between comparisons of samples, conditions, contrasts.

[18:01](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1081)
For the third branch that is based on DIA-NN, the data-independent acquisition branch, we made it fully parallelizable by a multi-step analysis. First you do an in silico library prediction and a pre-analysis for every sample, and then only after you do an empirical, data-dependent or data-based library generation, and a final analysis on the full experiment. It is also compatible with MSstats. This means you could have... the output will not be comparable in the quantities, but it will be comparable in the format compared to other branches of the workflow. As all the other branches, it can be converted into mzTab, which is, for smaller experiments, a human-readable tab-based format for the quantities and identifications of such experiments. You can use it immediately for upload to PRIDE, for example, or publication, which is usually recommended by the journals.

[19:36](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1176)
A bit more details on our general outputs. As I said, we have this mzTab for all the quantitative and identification-related information. The mzTab, in general, contains, on the right side here, metadata, a protein section, a peptide section, peptide spectrum match section, and for metabolomics also small molecule section. It's the community standard, so it's used by a lot of projects, and it's very helpful to have it for upload or journals.

[20:22](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1222)
Then, from our statistical post-processing, we can get heat maps or volcano plots for the comparisons between conditions that we can specify in the parameters of the workflow, for example. But we also have a full pMultiQC report, which is based on a plugin that we wrote for MultiQC, specifically for proteomics. It includes quality control heat map over all samples, but also detailed plots per sample and a detailed and searchable table of the results that is connected to an SQL backend. Those are some examples of our outputs. The first picture, you can see the experimental design that you have given and how it was interpreted by our tool. They can get very complex in proteomics experiments because you can also fractionate your data or your samples, and with the usual biological and technical replicates can get quite complex. In the lower part, you can see a heat map of some aggregated quality control metrics for specific samples and things like how many percent of contaminants were identified, the average peptide intensity, how many missed cleavages in your digestion we could find, what was the rate of identifications from the overall number of spectra, and so on.

[22:30](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1350)
Some more detailed information about specific samples. For example, the number of spectra on each level, like fragment spectra or service spectra, how many of them were identified by each of the search engines, how many were identified after consensusID, and so on. Lastly, one application that this workflow already had was reanalysis of a large part of PRIDE. We really sat down and were annotating with a large portion of PRIDE into the sample to data relationship format, which meant a lot of looking into papers, contacting authors, and so on. But which also means that you now, if you want to reanalyze something in a different way, you can just download the data from PRIDE or give URLs, which Nextflow, of course, handles to the FTP, and reanalyze it with different settings because the SDRF is already available.

[24:00](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1440)
We then reanalyzed each entry with our quantms. The good thing is we could analyze many of them because we made it very robust, our default settings, and also supporting a lot of different experiment types, as you have seen. Then in the end, we just combined and visualized the results, in this case per dataset or per tissue, because a lot of datasets are very specific for a certain tissue. We're currently writing a publication on that. That was one of the first applications, yes. I think that's it from our side. We're happy to answer all of your questions.

[24:56](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1496)
(host) Thank you so much. I'm just going to remove the spotlights. If there are any questions from the audience, you should now be able to unmute yourself and ask the question right away. Are there any questions from the audience?

(question) If not, I actually have a question, maybe a bit selfish. It's very nice to see that we have some pipelines at least that are not NGS based. I was wondering what made you choose Nextflow and nf-core for making this pipeline?

(answer) So the first thing was the incredible integration of all those large scale, high performance computing clusters and clouds that we have not seen in other workflow managers. Of course, a little bit bias because I knew some people from Nextflow. But I think it turned out to be the best choice in hindsight anyway. The nf-core team was very helpful in implementing all of this and the AWS tests were also super nice because as a university, we barely have any capability to test it on Amazon cloud or something that always costs. I think it gets a better reach also to industry by supporting clouds.

[27:02](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1622)
(question) Maybe in the same vein, did you find any problems that were specifically there because it is not NGS and because we're often very geared towards NGS?

(answer) Yes, of course. It's not big problems, but some of your templates, let's say, they have a lot of, not a lot... but what was it for example? I think you, or in the beginning you had a FastQC parameter that was always supposed to be there. We of course had to remove it. Now, whenever a template update comes, we have to remove it again and things like that. But yeah, minor things.

[27:48](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1668)
(host) Okay. Are there any more questions from anyone?

(question) I would have one. Hi. Great talk. I was wondering, you mentioned the small molecule MS experiments as a future possible application of the quantms. How far is this thought out or where does this stand?

(answer) Yeah. Implementation wise, we have a colleague that created such a workflow based on very similar tools that we already have. That means the OpenMS ones, but also some other tools like Sirius for small molecule database search. In the competitor language, Snakemake, but at least we see that it's a very feasible workflow that we have. It should be a rather simple translation of the workflow, but we also want to check with the existing MetaboIgniter workflow to see if we can combine them. We still have to check how compatible SDRF and the mzTab would be, so that we, for everything that we want to include into quantms we definitely want to start from an SDRF and have us output an mzTab or another future community standard file format. We think it should all be possible since there's also an mzTabM for metabolomics. Yasset, I think SDRF should have no problems at all to have some metabolomic specific annotations there.

(speaker 1) Yeah. I think we have a startup already to support metabolomics with SDRF. I mean, we have the first call around it, how to do it. As you said, I think this is a really important point. We have tried to put standard file formats in quantms as the starting point and the end of the workflow. For anyone who wants to join quantms this will be the case for other use cases like proteomics, like immunopatidomics or any other use case that wants to jump into mass spec quantitation in quantms. To start by one standard file format, something that the data out there is in that file format and should end up into another standard file format, which is this. In this case, it's mzTab, but it could be something in the future slightly different.

(comment) Great. Thank you both for the elaborations.

[31:18](https://www.youtube.com/watch?v=pBzelkgrPgQ&t=1878)
(host) Thank you very much. If we have any more questions from the audience? It doesn't seem so, then I would like to thank again, Julianus and Yasset, and of course, as usual, also the Chan Zuckerberg Initiative for funding our bytesize talks. Thank you very much everyone and see you hopefully.

</details>
