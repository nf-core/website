---
title: Hackathon - October 2021
subtitle: A virtual online hackathon to develop nf-core together
type: hackathon
start_date: "2021-10-27"
start_time: "10:00 CEST"
end_date: "2021-10-29"
end_time: "18:00 CEST"
location_name: Gather town and Slack.
---

# Introduction

> Our October hackathon is two weeks away! Registration is now closed, but keep an eye on your inbox if you've signed up, you will be hearing from us soon!

There will be a lot of firsts at this event. We are in the process of working out details with organisers arounds the globe to make it as geographically accessible as we can :earth_americas: :earth_asia: :earth_africa:. In addition, we will be experimenting with a group programming session on writing and contributing to nf-core modules for a limited number of participants. Our mentors, based around the world will be connecting with a small group of mentees who are also based in different countries!

This is a fair amount of work, and we want to ensure that you have a productive and fun event. So stay tuned.

Our primary focus for this hackathon will be the conversion of nf-core pipelines from DSL1 to DSL2. The main objectives of this hackathon will be adding nf-core modules and rewriting existing pipelines to the new Nextflow language format.

We have recorded bytesize talks in the past few months going over some of the details of tasks we will be tackling during the hackathon. Take a look if you  would like to learn more:

* [GitHub contribution basics](https://www.youtube.com/watch?v=gTEXDXWf4hE&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=4)
* [DSL module development](https://www.youtube.com/watch?v=ggGGhTMgyHI&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=5)
* [Adding modules to nf-core/modules](https://www.youtube.com/watch?v=Wc4A2tQ6WWY&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=7)
* [How to use modules in a pipeline](https://www.youtube.com/watch?v=tWvou0xj9wA&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=6)
* [Modules test data](https://www.youtube.com/watch?v=QXfAerydAT0&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=17)
* [Test modules](https://www.youtube.com/watch?v=pjhscKyWH74&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=18)
* [Development environments & workflows (Phil)](https://www.youtube.com/watch?v=XB96efweCLI&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=12)
* [Development environments & workflows (Maxime)](https://www.youtube.com/watch?v=OF55x-FT5WE&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=19)

## Hackathon groups

As the main objectives of this hackathon will be adding nf-core modules and rewriting existing pipelines to the new Nextflow language format, we have set up a few groups working on this task for different pipelines.

We will coordinate our work and the issues we are working on using a single GitHub [<i class="fab fa-github"></i> Project Board](https://github.com/orgs/nf-core/projects/22).

* **variant calling**
  * Working on adding modules for variant calling pipelines, mainly sarek, and raredisease.
  * Group lead: Maxime Garcia, Friederike Hanssen.
  * [<i class="fas fa-file-alt"></i> HackMD](https://hackmd.io/)
* **meta-omics**
  * Working on adding modules for meta-omics pipelines, mainly mag, bacass.
  * Group lead: Daniel Straub.
  * [<i class="fas fa-file-alt"></i> HackMD](https://hackmd.io/)
* **epigenetics**
  * Working on adding modules for epigenetics pipelines, mainly chip-seq, methylseq, atacseq, HiC.
  * Group lead: Harshil Patel.
  * [<i class="fas fa-file-alt"></i> HackMD](https://hackmd.io/)
* **ancient DNA (meta)genomics**
  * Working on adding modules for ancient DNA sequencing pipelines, mainly eager, coproid.
  * Group lead: James Fellows Yates.
  * [<i class="fas fa-file-alt"></i> HackMD](https://hackmd.io/)
* **single-cell**
  * Working on adding modules for single-cell pipelines, mainly scrnaseq, scflow.
  * Group lead: Gisela Gabernet.
  * [<i class="fas fa-file-alt"></i> HackMD](https://hackmd.io/)

## How we will work

We will be a lot of people working in parallel during this hackathon, so to stay organised we have a recommended workflow:

1. :speech_balloon: Chat with your group to get an overview of what is going on
2. <i class="fab fa-slack"></i> Join the relevant Slack channel to stay up to date and discuss with your project members
3. <i class="fab fa-github"></i> Find a task to work on using the GitHub Project Board
    * If you have something you want to do that's not there, please make an issue (e.g. in the nf-core/modules repository if you are adding a new module) and add it to the board
4. :raising_hand: Assign yourself to the issue that you're currently working on (preferably one issue at a time)
    * This is so that multiple people don't accidentally work on the same task
5. :fast_forward: When you're done, make a pull-request with your changes. Link it to the issue so that the issue closes when merged.
6. :page_facing_up: Describe your work on the HackMD document for the project and tell the group! :tada:
7. :recycle: Repeat!

> The HackMD document is the easiest to forget, but please add something even if you think what you did was small -
> we will use it in the group check-outs for each day, to stay tuned between time-zones and also in the reporting after the event so it's important for us :bow:

