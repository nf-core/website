---
title: 'Bytesize: Bactopia & using nf-core components in non-nf-core pipelines'
subtitle: Robert A Petit III - Wyoming Public Health Lab, USA
type: talk
start_date: '2022-07-19'
start_time: '13:00 CEST'
end_date: '2022-07-19'
end_time: '13:30 CEST'
youtube_embed: https://www.youtube.com/watch?v=egjgcmeJ0wQ
location_url:
  - https://www.youtube.com/watch?v=egjgcmeJ0wQ
  - https://doi.org/10.6084/m9.figshare.20338464.v1 (slides)
  - https://doi.org/10.6084/m9.figshare.20342808.v1 (video)
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: Bactopia & using nf-core components in non-nf-core pipelines

This week, Robert A. Petit III ([@rpetit3](https://github.com/rpetit3)) will introduce us to his pipeline [Bactopia](https://github.com/bactopia/bactopia). A first for bytesize, bactopia is not an nf-core pipeline but is rather [influenced by nf-core](https://bactopia.github.io/acknowledgements/#nf-core).

Robert will show us the pipeline and explain how he is using nf-core components (such as [modules](https://nf-co.re/modules)) in this pipeline, despite it being outside the main nf-core organisation.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=1)
(host) Hello everyone my name is Franziska Bonath, I'm the today's host for the bytesize talk and with me is Robert Petit. He is from the Wyoming public health laboratory and is going to talk about bactopia which is not part of nf-core but it's using nf-core components and nf-core pipelines, so off to you!

[0:28](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=28)

(host) I'm so sorry I did mute you by mistake. You can unmute yourself. I'm sure...

There we go!

(host) I'm so sorry!

No don't worry about it!

[0:50](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=50)
Thank you for having me. I'm super excited to be here today. I'll just be introducing you to bactopia. I'll give you a few use cases of how I'm using nf-core components in bactopia. But first a little motivation on how bactopia came up. Over the last 10 years we saw a a quite nice growth in the publicly available bacterial genomes. From ENA/SRA/DDBJ we went from about 7500 in 2010 to about 1.5 million. While that pales in comparison to COVID, which I think is 12 million plus now, for bacteria that is quite a lot. Over these last 10 years we've also seen the rise of containers and package managers, such as docker, singularity and bioconda. And then workflow managers and the curators behind those, such as nextflow and nf-core. At 2010 maybe we couldn't use all the data but it really starts in 2022, it really makes you wonder: can we make use of all these publicly developed genomes? In 2010 I remember, passing tar balls with binaries and then emailing the sequencing instrument groups say "can I get a binary to your assembler?" and stuff like that. There was no real... It was just "make install" and hope it installs.

[2:28](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=148)
Now in 2022 I think we have the tools and we have the talent to start gluing all this stuff together and start using these data in our own analyses. Once you know why would we want to use this data. A good example is, if you have a small outbreak at your local hospital, or a foodborne illness that comes out of some carnivor and you want to compare your genomes to what's already been sequenced. That's a nice use case of making use of those 1.5 million genomes that are available. To address this, Tim Read, who's at Emory University and was my master's and phd mentor, him and I develop bactopia, which is a nextflow DSL2 pipeline for the complete analysis of bacterial genomes. Because it's written in nextflow you can go from a single genome to tens of thousands of genomes with a simple parameter change.

[3:33](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=213)
We were able to use bactopia to process six to seven thousand genomes in just five days using AWS. A lot of that was being able to prototype on a laptop and then switch to our AWS profile and "boom!", we're off and going. Kudos to the nextflow for most of that. In bactopia we try to include as many nf-core practices as we can, to ensure things like reproducibility and audit logs and all that. Because it's nextflow it's extremely portable and you can go laptop, HPC at your university or something between all the cloud platforms within just a few parameter changes.

[4:28](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=268)
A few highlights about bactopia. It supports Illumina and nanopore reads. These can be from your local machine or from publicly available databases, such as SRA or ENA. It includes more than 140 bioinformatic tools. There are 45 bactopia tools, which are completely separate workflows, which I'll get into shortly. It's been extensively tested with more than 100 tests, tested more than 10.000 output files. It's easily installed through bioconda, docker, or singularity and I've gone through great efforts to make sure it's well documented.

[5:09](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=309)
Some design principles behind bactopia. One bactopia requires all tools that are included in it to be available from bioconda. The main reasons is because it's 2022. People shouldn't have to figure out how to install a tool now. It should be able to either use a container or some sort of conda. It should be an easy, simple process and because bioconda has the downstream containerization, so every recipe gets a docker container through bio-containers and a singularity image through galaxy project, we have all those tools necessary to start using these immediately. I also require all modules and bactopia tools to also be available from nf-core modules. If it's not there we add it. Bactopia should be easy to install and adaptable to the user's needs. Converting to DSL2 has made this much, much easier.

[6:14](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=374)
There are three sides to bactopia. You can think of these as checkpoints between the three. There's bactopia helpers, bactopia and bactopia tools. The bactopia helpers help you get started using bactopia. These are your pre-analysis steps or some commands to post analysis, get information. One is the bactopia citations, which will print out citations for all the tools used by bactopia. The bactopia datasets command allows you to go and download publicly available data sets, that can supplement your analysis. These include things like RefSeq and GenBank sketches, as well as pub MLST schemas and many more. The bactopia download command will pre-build conda environments for you, pull docker containers, or download singularity images as a pre-step so that you're not doing that while you're starting to process the nextflow. The bactopia prepare will create a file of file name similar to the sample sheet that you see in many of the nf-core pipelines. This allows you to really process as many genomes as you want. The bactopia search, one of my favorite, it takes a query, queries ENA's API, then returns a list of experiment accessions that you can then feed to bactopia to download and start processing.

[7:49](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=469)
The main bactopia pipeline includes all the standard steps in a bacterial genome analysis: Gather samples, QC the reads, the simple genomes. You can sketch your genomes and then query against RefSeq or GenBank, call SNPs, ... All the standard things that you would expect in a bacterial genomics pipeline. It allows Illumina or nanopore reads, SRA accessions, NCBI accessions or local assemblies, if for some reason that's all you have. There are all also some jump off, where basically the sample will stop being processed, if there's things such as poor quality, something that's gonna likely cause downstream failure. Bactopia will do its best to catch those, so that way it doesn't stop the whole pipeline. Once everything's processed you get it in this nice standard directory structure. It's this directory structure that bactopia tools take advantage of. Bactopia tools are essentially more workflows for more science. By looking at that standard directory structure, you can run a bactopia tool, which can include a single tool like Kleborate or TB Profiler, and then it'll go and find the files that it needs and run everything for you. You can connect multiple modules together for something like a pangenome type analysis, where you're running PIRATE and creating a core genome phylogeny. The bactopia tools, because of the directory structure, will find all those files that you need. There are currently more than 45 different bactopia tools. Because it is DSL2 I've been able to framework this and make it a streamlined process.

[9:46](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=586)
In just a few steps you can go from raw data to investigating results: (1) Sequence your genomes, (2) install bactopia through conda, docker, singularity, (3) if you want to include public data you can use bactopia search. (4) If you want to include publicly available data sets, which I always recommend to supplement your analysis, use the bactopia datasets command, and then (5) you can create a file of filenames to process thousands of genomes, if you want to using bacteria prepare. You use (6) the bactopia command to process all your samples independently and then (7) further analyze these with bacteria tools. (8) By the end of it you just have a bunch of output files that you get to sift through and figure out: can we answer our question, that we hopefully asked before sequencing these genomes.

[10:43](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=643)
Most of this has been made easier and more achievable in bactopia by adopting nf-core components. If you're on the outside wondering, should I make an nf-core pipeline or should I just keep doing what I'm doing and start adopting some of their practices, or should I just go do my own thing, I think you're the target audience here. Over these next few slides you can get an idea of how I am making use of numerous nf-core components, without actually being an nf-core pipeline. Honestly, I don't think it's so much about the nf-core practices and components and more about the people behind nf-core. You jump on the slack group, and you got a question, and there's many people that are willing to help out. They've probably seen it, especially many of the error messages that you come about in bioinformatics. I think at minimum you should hop onto the slack group and just start participating and get an idea of all the things happening with nexflow.

[12:02](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=722)
Here are a few ways I'm making use of nf-core components in bactopia. First the nf-core library, which is that lib folder in all the nf-core pipelines. Bactopia has 45 different workflows that you can execute from a single entry point. So there's a parameter that says "I want to run the pangenome bactopia tool" or "I want to run the bactopia main workflow". Those all come in through the same main .nf file. To achieve this I adapted the nf-core library, because 1, it handles all the argument parts, and it has super nice outputs. It does auditing and you can set it up to send emails and all that. Also by using that, you set yourself up to be compatible with nextflow tower, which is quite nice.
But I wanted to be able to programmatically import config and json files. On the bactopia side, I have a dynamic import, that looks at a workflow config and determines based on that, which files it needs to import. That way I can run 45+ different workflows from the same entry point, which is quite nice, because previously it would have been 45 different main .workflows that I was maintaining.

[13:32](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=812)
Nf-core modules in bactopia. When I converted DSL2, it was suggested to me that I should consider making use of nf-core modules. I had previously participated in some of the hackathons and was quite fond of nf-core modules, so it was super easy to say, okay, if I'm going to include a bactopia tool, it should also be on nf-core modules.

[13:56](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=836)
On the bactopia side, I do some slight modifications. These slides will be available later. There's links to compare the two, there's many links in in these slides. Some of those modifications are mostly just adapting to use pre-built conda environments. Just the way I import and export files. I've also adopted a similar pytest framework for bactopia, that is implemented in nf-core modules. This allows me to test every step in bactopia and bactopia tools. This has saved my butt quite a bit, when it comes to submitting a new release. Typically it's the conda side where something has changed with the package solver. This also allows me to use self-hosted github action runner. Those modules, like gtdb, which use large databases, are being tested with a real database on my self-hosted github action runner. There's that side effect that we're validating indirectly the nf-core module.

[15:05](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=905)
Finally, I used the the meta .yaml template for documentation. When I first saw that meta .yaml I was "oh, that would be nice to just build documentation from". I add stuff like citations, some Markdown tables and output trees. The YAMLs are then used to build the documentation using Jinja2 templates, MkDocs material and GitHub actions. This has really saved me a lot of time, by allowing me to write the documentation while building bactopia.

[15:38](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=938)
What's next for bactopia? I am always waiting to see what's next for nf-core and in the background saying "hey, should I use this or not?". I'm starting to look into multi-qc modules, because bactopia needs some sort of report generation. There will always be more nf-core modules that I'm submitting, because there's always more bactopia tools I want to submit or implement. I really have my eyes on the that issue on nextflow about the future of the config files, because the way I use config files. That could have some some downstream effects on bactopia. I'm interested in making a custom workflow for surveillance, here at bactopia. The more I use rich-click I just want to rich-click everything. Expect an enhanced cli here soon.
Don't hesitate to reach out if you think I can help you get started on your non nf-core pipeline and using nf-core modules.

[16:37](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=997)
Thank you and I will take any questions folks have.
(host) Thank you very much. I have enabled now for everyone to unmute themselves, if they have questions. Otherwise we can start with the one question that is already in the chat, which is from Olaitan. He is asking: "it seems like bactopia is not in the APT repository, could you work on including it?"
(answer) I don't know much about including tools in the APT repository. It is available from bioconda, so you can "conda install" bactopia. I think there would be many components of bactopia, that aren't in the APT repository, so I don't know how that would work. It would be something you would have to add all the dependencies to the APT repository, and I think the time required for that and the learning, I don't have the bandwidth for at the moment. Definitely consider using the bioconda install and then from there you can use conda, docker, singularity...

[17:55](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=1075)
(host) okay thank you. Are there any more questions from the audience?
(question) I do have a follow-up question. It's not really based on the APT repository thing, it's the workflow. In terms of the different steps, that bactopia does. I didn't get what Robert meant by the final step that he talked about. Something about the analysis. I was wondering what are the specific things, like what are you measuring at the end of the day with bactopia? Specifically in terms of the omics analysis.
(answer) It's going to include pretty much all your standard bacterial genomics. You're going to QC the reads, how well did you sequence your sample, what's the average read length, all that fun stuff. Then it's going to characterize your sample. what MLST schema? Does it have certain antibiotic resistance genes that you may be interested in? Does it have SNPs and INDELs against a reference genome that you selected? How does it compare to publicly available genomes, does it look like what you expected? If you thought you sequenced the Staphylococcus aureus, and you know it came up as looking more like Enterococcus, that's something. Those are the type of analysis results. On the bactopia documentation there's an overview of the workflow at each step and then output overviews on all the output files that you get for both bactopia and all the bactopia tools. Those output give you a description of what's in the each of the files.

[19:52](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=1192)
(question) Okay thanks! My final question is, I saw you integrated this with Illumina reads and I believe Nanopore reads. What happened to PacBio?
(answer) Honestly I just haven't been exposed to PacBio data much so far in my analyses and my studies. I think if I start using PacBio, then PacBio data will come in, otherwise I think I would need support from the community to add that type. Just because I don't have the opportunity to use it on a daily basis like I do Illumina and Nanopore. We need someone else to help out there.
(from audience) cool, thanks! Great job! thank you!
(host) thank you very much, are there any other questions?

[20:55](https://www.youtube.com/watch?v=egjgcmeJ0wQ&t=1255)
(host) I don't see anything pop up so, I would like to thank you again, Robert. Everyone else, there is always the chance to ask more questions, if they come later, at the bytesize channel on Slack. I guess you can also contact Robert directly and this video will also be uploaded to Youtube. I would like to thank, apart from Robert of course, the Chan Zuckerberg Initiative, who is funding these talks, and thank you everyone for joining in. Thank you.

</details>
