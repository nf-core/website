---
title: new modules to infer protein activity
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C08K6ACBUHH
image: "/assets/images/events/2025/hackathon-march/protein-activity.jpg"
image_alt: "Buzz looking at the protein activity interactions to the infinity and beyond"
leaders:
  arangooj:
    name: Juan E. Arango Ossa
    slack: "https://nfcore.slack.com/team/U04LMD2LPB6"
  domenicd:
    name: Dylan Domenico
    slack: "https://nfcore.slack.com/team/U08JZQY91JR"
  liosisk:
    name: Konstantinos Liossis
    slack: "https://nfcore.slack.com/team/U08J7L7E0T1"
---

Tools like [VIPER](https://static1.squarespace.com/static/5697c2e5e0327ca6778bc453/t/56f40934f8baf3727f8e7e78/1458833718573/Viper.pdf) and [ARACne](https://github.com/califano-lab/ARACNe3) from [the Califano Lab](https://califano.c2b2.columbia.edu) are useful in understanding the protein-activity interactions, on an individual sample basis, from gene expression profile data:

"Under the premise that cellular phenotypes not only result from alterations in the genomic code but also depend on the influences of multiscale networks of molecular interactions that regulate gene expression, protein abundance, epigenetic state, and signaling activity. Therefore, understanding how the code embedded in the genome eventually generates cellular phenotypes requires the development of methods for integrating data describing these various levels of activity." [Source](https://califano.c2b2.columbia.edu)

### Goal

Create modules and a subworkflow for the computational inference of protein activity network and its interactions, useful for precision medicine to understand how drugs are involved in protein pathways.
