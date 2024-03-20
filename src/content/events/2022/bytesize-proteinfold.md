---
title: 'Bytesize: nf-core/proteinfold'
subtitle: Athanasios Baltzis - Centre for Genomic Regulation, Barcelona
type: talk
startDate: '2022-09-20'
startTime: '13:00+02:00'
endDate: '2022-09-20'
endTime: '13:30+02:00'
embedAt: 'proteinfold'
youtubeEmbed: https://www.youtube.com/watch?v=441P_GzrI-o
locations:
  - name: Online
    links:
      - https://www.youtube.com/watch?v=441P_GzrI-o
      - https://doi.org/10.6084/m9.figshare.21206555.v1
---

This week, Athanasios Baltzis ([@athbaltzis](https://github.com/athbaltzis)) will talk about the newest developments in the nf-core/proteinfold pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=441P_GzrI-o&t=1)
Hi Maxine here. Today I'd like to welcome Athanasios Baltzis from the comparative bioinformatics group at the center for genomics regulation in Barcelona. He's going to talk about proteinfold, which is pipeline I'm really looking forward to know more. It is before the first release, which is from what I hear, coming soon. Before we actually start I'd like to thank thank the Chan Zuckerberg Initiative for helping us out, and the listeners. You will be able to admit yourself at the end of the talk for questions.

[0:42](https://www.youtube.com/watch?v=441P_GzrI-o&t=42)
thank you Maxime for the nice introduction. I'm very happy today that I will present the nf-core proteinfold pipeline, best practice bioinformatics pipeline for protein structure prediction.

[0:58](https://www.youtube.com/watch?v=441P_GzrI-o&t=58)
Let me introduce first myself.I'm a PhD fellow in Cedric Notredame's Lab at the Center for genomic regulation in Barcelona and my thesis is about applications of protein structure modeling on multiple sequence alignment and phylogenetic reconstruction. This is where it connects to protein structure prediction methods because I am very interested in all these tools. I use them in my daily routine. Let's first start with the brief introduction to protein structure prediction. As you may know there is cast data on experimental protein structures, mainly due to technical difficulties with the already existing techniques. It is a long-standing question for the community how can we start from 1D amino acid sequence to go to a 3D structure and gain more insight into the function of the proteins, this is the the so-called protein folding problem. For this reason many techniques have been developed during the the last mid-century and can be categorized into two main categories: the template based methods ,for example homology modeling or fault recognition which rely on already existing experimental structures that are used as templates to fold the query sequence. And on the other hand we have also template free ab initio methods. We have a lot of categories, for example molecular Dynamic simulations, where we try to use physics law to find the confirmation with the lowest dips free energy. There are fragment-based approaches such as Rosetta and lately pairwise special restraint based approaches where you use techniques to predict the inter-residue contacts or the distances between the query protein and then use them as restraints in simulations, in order to get the final predictive model.

[03:49](https://www.youtube.com/watch?v=441P_GzrI-o&t=229)
In the last year, AlphaFold2 achieved a major breakthrough in this field and it's now able to predict protein structures from sequence with an unprecedented accuracy, near experimental accuracy, I would say. This is mainly based on the incorporation of deep learning frameworks in the field. For example here in the figure below, you can see a brief representation of the AlphaFold workflow where we start with an input sequence. You search genetic databases in order to form an MSA from a homologous sequence and convert it into a tensor. On the other hand if you search for structural templates in order to populate this pairing Matrix, which actually represent the interest interactions of the input sequence. AlphaFold consists of two main neural network blocks: the Evo former, where it gets ensembled the MSA representation and the pair representation. The MSA representation tries to populate and optimize the pair representation Matrix. Afterwards we have the structure module, where we convert these two type of tensors into a tensor that contains the translations and rotation of the model and during the learning process this is optimized. Finally it gives us a final 3D structure. Of course for better performance and accuracy this happens three times, there are three recycling steps. After the release of AlphaFold there were several other tools with similar, or even better, accuracy and performance than AlphaFold. But the problem with this software is, that they have a lot of dependencies, and mainly we refer to the genetic databases you need at this step here, in order to build the input multiple sequence alignment.

[06:44](https://www.youtube.com/watch?v=441P_GzrI-o&t=404)
As many labs and researchers in the community try to use Alpha fold in a large scale - as we did in our msaaf2 nextflow pipeline - we were interested into develop a pipeline that can take care of all these dependencies: the databases, the AlphaFold parameters,... in order to be able to get fast and as reproducible as possible predictive models. After the release of our pipeline here many researchers got in contact with us, from the academic or the private sector, that were interested for a scalable AlphaFold pipeline, that deals with this problem with the dependencies. For this reason we got in contact with nf-core and Seqera labs and we collaborated in order to develop such a scalable protein structure prediction pipeline.

[08:03](https://www.youtube.com/watch?v=441P_GzrI-o&t=483)
Here we have an overview of what we already have at the moment. As you can see we have four modes, mainly based on two sub workflows: The AlphaFold2 and the CollabFold. Let's go through this overview step-by-step. We start with the input sample set which is quite similar with the the input that is already used in the majority of nf-core pipelines. It's a bit different in the sense that here we have a comma separated file with two columns. The first column is the sequence header and the second column is the path to the .fasta file. For monomer predictions it is recommended to use multiple entries. For its monomer sequence you want to predict and here you have an example of a fasta file. For multiple predictions it is recommended to use one entry with a corresponding to a multi-fasta file. For example here you have this multimer and you have a look here at the multifasta file containing as many entries as the sub units that you want to predict for this multi-mer.

[09:38](https://www.youtube.com/watch?v=441P_GzrI-o&t=578)
Once the pipeline checks for the validity of the input, you have two options, two sub-workflows. The first one is the AlphaFold2. The sub workflow where it first passes through the prepare AlphaFold2 sub-workflow, which checks if you have provided the af2db parameter, which specifies the path where the pipeline can find all the databases and the parameters that AlphaFold will use, if you have downloaded them. Otherwise it downloads those themselves to the required databases and model parameters. I would like to point out at this point that this is quite computationally expensive since it has to download around 2.2 terabytes.

[10:42](https://www.youtube.com/watch?v=441P_GzrI-o&t=642)
You can use AlphaFold in two modes. The first one is the default one where you just feed the input CSV to the AlphaFold and it computes the input multiple sequence alignment and does the model inference in the same process. But for the sake of computational cost, we also implemented another mode of AlphaFold, we call it AlphaFold split, where it gets the input CSV and the .fasta and produces the input MSA in a separate process than the model inference. If you think about it, this is quite convenient for example in Cloud infrastructures because this step, the af2 prediction step, requires GPU. If you run these two steps in the default mode on GPUs, it's much more computationally expensive, it costs much more than af2 split. You can specify these two modes using the standard, the af2 parameter here, TRUE for the default AlphaFold and FALSE for the AlphaFold split.

[12:18](https://www.youtube.com/watch?v=441P_GzrI-o&t=738)
The second sub-workflow is about COLABFOLD. We have more or less a similar strategy. We have the prepare COLABFOLD sub-workflow, where you can specify if you have already downloaded the databases and the required parameters of the model. You can specify the path using the COLABFOLD DB parameter, otherwise it downloads automatically the required databases and model parameters and here again it requires a lot of storage - around 1.8 terabytes. We have two modes in the COLABFOLD as well. The default mode is the COLABFOLD web server, where you depend on a web server that can run the database search and MSA creation. By default this web server is the one provided by a mm6, team but using the parameter host URL you can specify the URL to your custom web server, f you have set it up. In order to to specify this mode - the COLABFOLD web server - you just have to use the mode parameter. The second one is the COLAB-for-local mode and it's quite similar to the AlphaFold split mode we have seen in the last slide, in the sense that if you first have a process to compute, the input multiple sequence alignment (we're using the mm6)and then you have a separate process for the model inference and the protein structure prediction. You can specify it using the mode parameter of the pipeline.

[14:20](https://www.youtube.com/watch?v=441P_GzrI-o&t=860)
Let's have a look at the some more advanced parameters. For example the --use_gpus parameter when available, because as I explained before it's a much more computationally expensive to use only CPUs, especially in the prediction steps. But you should also pay attention to the configuration profile we are using in combination with the GPU. To define the corresponding prediction process to the GPU or machines you have in your infrastructure. For example we have in the GitHub repository of the pipeline a CRG institutional profile that we are using at the moment, so we can have a look. With the --outdir parameter you can specify the output directory. This applies for all of the nf-core Pipelines.

[15:38](https://www.youtube.com/watch?v=441P_GzrI-o&t=938)
Some specific AlphaFold2 options. The --full_dbs parameter where you can select if you want to use the full databases for the first part of the sequence search and the creation of the MSA, or you can use a reduced version of the databases, which means that the pipeline will run faster, but with a bit of a trade-off in accuracy of the produced model. The --model_preset parameter, where you have to specify what type of prediction you want to do and which model to use. For example we have a three monomer model, the default is this one. The other two are actually provided by the AlphaFold2 team for reproducibility purposes. This one was used in the Casp competition casp14 competition for example, or the multimer model for multiple predictions.

[16:52](https://www.youtube.com/watch?v=441P_GzrI-o&t=1012)
Regarding the COLABFOLD specific options. You can specify the model type, the AlphaFold PDM, which is the default for monomers and to a multimer models the the most improved version is the default the version 2. You can also specify if you want to use pdb structure templates or not in the first step, where you populate the pair representation Matrix. You can find also some more specific and detailed description on the parameters available at the moment in the corresponding web page of the nf-core proteinfold pipeline.

[17:43](https://www.youtube.com/watch?v=441P_GzrI-o&t=1063)
Regarding the output. Here you can see the tree structure of the produced output. If you use AlphaFold you have an AlphaFold directory and one more subdirectory with the sequence name you have provided in your CSV file that contains the computed MSA. It further contains the unrelaxed and relaxed structures the rank structures the raw model output, some metadata and of course timings. The first ranked model that probably is what you want to use in your research and that contains as well the pldt scores, which is the confidence metric used by AlphaFold per residue. Another directory contains symbolic links to the downloaded databases and parameter files. The same applies for COLABFOLD, where you have an output directory depending on the mode you have selected: COLABFOLD web server or COLABFOLD local, that contains all this information we have explained for AlphaFold, and the symbolic links to the downloaded databases. Of course as in all the nf-core pipelines there is a directory with the pipeline in for execution Trace files and so on.

[19:28](https://www.youtube.com/watch?v=441P_GzrI-o&t=1168)
What are the next steps? We are now at this point that we have to set up and run the AWS full tests, in order to create the first release of the pipeline. In future releases we are planning to add more open source proteins to actual prediction tools, such as open fold, or even a newer generation of prediction tools, such as esm fault or Omega fault, which use protein language models and are for this reason about an order of magnitude faster than AlphaFold or COLABFOLD, without losing accuracy. In fact they have the the same levels, or even better levels, of accuracy. We're also interested in incorporating more advanced software for protein-protein interaction such as FoldDock, because there are plenty of researcher interested in predicting a advanced multimers, and moreover, to add to solve bug fixes and add more optimization upon request, we are very open to contributions and ideas, in order to improve even more the pipeline and adapt it to the needs of the community. So please do not hesitate to contact us and propose or contribute to the already existing Repository.

[21:27](https://www.youtube.com/watch?v=441P_GzrI-o&t=1287)
At this point, before finishing, I would like to thank my colleagues. First from the Notredame's lab, Jose Espinosa and Luisa Santus, that are contributing to this pipeline, as well Seqera labs and especially Harshil Patel, for all the guidance and the help during the implementation process. Also I would like to thank the collaborators from the Interline Therapeutics, especially Norma Gudager and Walid Usman, for testing the pipeline in the cloud.

[22:07](https://www.youtube.com/watch?v=441P_GzrI-o&t=1327)
Thank you for your attention and I would be very happy to answer if you have any question and that's it, thank you.

(Maxime) Thank you very much for the amazing talk. I will allow everyone to unmute themselves if anyone has any question. Please, let's go.

[22:34](https://www.youtube.com/watch?v=441P_GzrI-o&t=1354)
(Question) Otherwise I have one question. So at the moment you only have AlphaFold2. And you are planning to add more tools, but not in this first release but in the comming one, right? I assume that the main issue with having more tools is that, it's a lot of databases that you need as an input.
(Answer) Exactly that's true because each tool uses its own databases, so you need a lot of storage to be able to test everything or even to compare between tools.

[23:16](https://www.youtube.com/watch?v=441P_GzrI-o&t=1396)
(Question) May I ask along this line. So you basically retrain the model every time you run the pipeline, or at least like every time an institution retrains their model from scratch, or do you use pre-trained models.
(Answer) We use pre-trained models. We just download the already provided models by AlphaFold
(Question) And it still takes these huge databases?
(Answer) yes, because this is separate from the training process. These databases are needed in order to create the input multiple sequence alignments, to actually have this or all these bunch of homologous sequences, in order for the model to be able to find all the correlations, the interesting new correlations and form the final model.

[24:15](https://www.youtube.com/watch?v=441P_GzrI-o&t=1455)
(Maxime) I think we are good with the number of questions. Thanks again, that was an amazing talk now I'm super happy to have learned more about it. I'm really hoping like to see this release coming.

</details>
