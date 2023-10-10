---
title: 'Bytesize 21: nf-core/mhcquant'
subtitle: Leon Bichmann - University of Tuebingen, Germany
type: talk
start_date: '2021-09-28'
start_time: '13:00+02:00'
end_date: '2021-09-28'
end_time: '13:30+02:00'
embed_at: 'mhcquant'
youtube_embed: https://youtu.be/NCKkSssE_4w
location_url:
  - https://youtu.be/NCKkSssE_4w
  - https://www.bilibili.com/video/BV1Lh411J732/
  - https://doi.org/10.6084/m9.figshare.16750381.v1
---

# nf-core/bytesize

Join us for a special pipeline-focussed episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 21: nf-core/mhcquant

This week, Leon Bichmann ([@Leon-Bichmann](https://github.com/Leon-Bichmann/)) will tell us all about the nf-core/mhcquant pipeline.

nfcore/mhcquant is a bioinformatics analysis pipeline used for quantitative processing of data dependent (DDA) peptidomics data. It was specifically designed to analyse immunopeptidomics data, which deals with the analysis of affinity purified, unspecifically cleaved peptides that have recently been discussed intensively in the context of cancer vaccines.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:41](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=41) Thanks for inviting me to present the nf-core/mhcquant pipeline that I developed during my PhD. nf-core/mhcquant is an automated pipeline to analyse mass spectrometry data for the discovery of epitopes that can be used for the design of vaccines.

[1:20](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=80) I’ve structured this talk to consist of three parts. I will first provide a short context by recapitulating cancer immunotherapy and mass spectrometry, then I will go into more details of the pipeline itself, and finally I will provide some future perspectives.

[1:49](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=109) This is a very basic image of T-cell-based adaptive immunity, which is just one branch of the immune system, but that was the part we were interested in. T-cells check all the cells in our body by checking one cell-surface protein complex called MHC (Major Histocompatibility Complex). This complex presents small peptides on its surface that represent fragments of all the protein content of a given cell. So if the cell is unhealthy, say due to a viral infection or a cancerous protein, then the T-cells can recognise this from the peptide epitopes and commence a cytotoxic activity that can kill the malignant cell.

[3:14](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=194) So this process is being exploited for cancer immunotherapy by comparatively analysing tumour and normal tissues using sequencing tools and mass spectrometry to investigate the MHC epitopes presented by the tumour tissue. Those then serve as candidates for vaccines that would stimulate the patient’s immune response against the tumour. The nf-core/mhcquant pipeline focuses on the identification of these epitopes.

[4:05](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=245) To give you a bit more of a feel for the kind of data we deal with, one takes tissue samples in the lab, homogenises them, purifies the MHC complexes and upon elution of the peptides from the complex, one can spike these into the MS instrument and measure this.

[4:37](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=277) Mass Spectrometers are high-throughput instruments so one can automatically sample from a box containing vials to obtain hundreds to thousands of runs in short timeframes.

[5:01](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=301) We therefore need high-throughput computational analysis pipelines to process the data that results from these runs.

[5:14](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=314) So this is where the nf-core/mhcquant workflow came into play and I’ll give you an overview of what is happening inside the pipeline.

[5:23](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=323) If you have a bit of an overview of the architecture of all nf-core pipelines, you know that there are processes carried out by software libraries at the centre - in this case the [OpenMS software library](https://www.openms.de/) for computational mass spectrometry. We combined this with third party tools and scripted things that weren’t available using python and R. This was all containerised using Docker and Singularity or other methods provided by the nf-core template and was run using the Nextflow workflow system in a highly reproducible manner.

[6:30](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=391) So here you see a rough sketch of the pipeline with its 35 different steps that are all interconnected in different ways.

[6:50](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=410) Let us focus on the five most important ones, and I will guide you through them step-by-step.

[6:55](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=415) First a protein database search is carried out using the search engine [Comet](http://comet-ms.sourceforge.net/).

[7:09](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=429) We used this well-established search engine because it is a very simple and fast-scoring method, and it has no tryptic bias for unspecifically cleaved peptides, which is what the MHC peptides are. In a benchmark comparing a variety of different tools, it appears that [Comet](http://comet-ms.sourceforge.net/) finds most peptides on average. The only other one that compares is [Peaks](https://www.bioinfor.com/peaks-studio/), however, that is a licenced tool and cannot therefore be used in this case.

[8:04](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=484) We also verified these additionally identified peptides comparing their retention time properties, and found that peptides identified by the Comet MHCquant workflow nicely correlated when compared to random decoy peptides that scattered all over the place.

[8:41](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=521) So as a next step, a classical thing in proteomics is the false discovery rate (FDR) and we use an advanced method called Percolater. For details, please consult [Käll et al., 2007](doi: 10.1038/nmeth1113.). In contrast to the classical approach of simply computing the FDR by looking at univariate target decoy distribution, a multivariate distribution is achieved by an iterative machine learning approach where the spectrum matches are compared with a variety of different scores. This results in a better discrimination between the target and the decoy.

[9:44](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=584) So in the next step, the retention time for all the peptides is corrected because the retention time can be slightly different across different measurements. This has been corrected in the pipeline using a specific tool called MapAlignerIdentification within the [OpenMS](https://www.openms.de/) tool. This aligns all the peptides to one central reference.

[10:30](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=630) So finally, we get to the part that explains were MHCquant gets its second name. Every peptide is associated with a quantity. This is carried out using a targeted chromatogram extraction approach - each peptide identification is located not only on the MS2 level but also at the MS1 level. The corresponding chromatograms are integrated and the sum of the signal intensity area under the curve represents the quantity that can be compared.

[11:11](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=671) Again, we went into the lab and also verified this so looking at the signal intensities of 57 spiked-in peptides that were diluted in a linear series, we observed a linear decay in signal intensity. So again, we validated that the quantification works reliably well.

[11:34](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=694) Finally, an MHC affinity prediction is carried out, and here we applied two open-source neural network architectures called MHCFlurry and MHCNuggets. We were very happy with how they work.

[12:17](https://youtu.be/NCKkSssE_4w?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=737) So with that, I’ve come to the end of my presentation but I’d like to give you some future perspectives. The pipeline has been ported to DSL2 by [Marissa Dubbelaar](https://github.com/marissaDubbelaar). The identification rates can be improved by intensity prediction, and since mass spectrometry instruments are constantly being improved, it would be nice to include ion mobility data and data-independent acquisition-based methods. Thank you for your attention. You can contact me on Slack if you have any questions.

</details>
