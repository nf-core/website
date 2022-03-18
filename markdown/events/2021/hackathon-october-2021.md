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

Our October hackathon is one week away! Registration is now closed, but preliminary information should be in your inbox if you've signed up. Final details will be coming soon.

There will be a lot of firsts at this event. We are in the process of working out details to make it as geographically accessible as we can. In addition, we will be experimenting with a group programming session on writing and contributing to nf-core modules for a limited number of participants. Our mentors, based around the world will be connecting with a small group of mentees who are also based in different countries! As always, there will be a few social activities throughout the event as well as on the Thursday evening (CEST).

# Location

This edition of the nf-core hackathon will be hosted on [gather.town](https://www.gather.town/) so please [familiarise yourself](https://support.gather.town/help/movement-and-basics) with the platform in advance. The link to the space will be sent a couple of days before the event.

Please also join the [nf-core Slack workspace](https://nf-co.re/join), and find us on the #hackathon-oct2021-public channel

# Prerequisites

Prior the hackathon, make sure you're signed up/joined/have installed the following resources necessary for participating in the event:

- Check you agree with the [Code of Conduct](https://nf-co.re/code_of_conduct) of the event.
- If you haven’t already, set-up a GitHub account and join the nf-core GitHub organisation.
- Join the [nf-core slack](https://nf-co.re/join) and the channel #hackathon-oct2021-public
- Have installed on your computer:
  - [Nextflow](https://nextflow.io/)
  - [nf-core/tools](https://nf-co.re/tools)
  - Docker/Singularity/Conda: [Google is your friend]
- Familiarise yourself with the documentation on the nf-core website for nf-core modules:
  - [https://nf-co.re/developers/modules](https://nf-co.re/developers/modules)
  - [https://www.nextflow.io/docs/latest/dsl2.html#modules](https://www.nextflow.io/docs/latest/dsl2.html#modules)
  - Relevant nf-core/bytesize talks are also listed below
- Have a peek at the [GitHub Projects board](https://github.com/orgs/nf-core/projects/20) for the hackathon
- If you have been accepted as a mentee for the pair programming session on day 1 your mentor should be in touch via DM on Slack. If you haven’t been contacted please let us know on the #hackathon-oct2021-public Slack channel.

If you have any problems with any of these just ask on the slack channel or email [outreach@nf-co.re](mailto:outreach@nf-co.re)

# Aims and Objectives

Our primary focus for this hackathon will be the conversion of nf-core pipelines from DSL1 to DSL2. The main objectives of this hackathon will be adding nf-core modules and rewriting existing pipelines to the new Nextflow language format. To get started we added two tutorials:

- Current state of the DSL2 template:

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/0xjc7PkF1Bc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>


- Tutorial: Adding a local module to nf-core/modules:

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>

<div class="ratio ratio-16x9">
    <iframe src="https://widgets.figshare.com/articles/16825369/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>

We have recorded bytesize talks in the past few months going over some of the details of tasks we will be tackling during the hackathon. Take a look if you  would like to learn more:

- [GitHub contribution basics](https://www.youtube.com/watch?v=gTEXDXWf4hE&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=4)
- [DSL module development](https://www.youtube.com/watch?v=ggGGhTMgyHI&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=5)
- [Adding modules to nf-core/modules](https://www.youtube.com/watch?v=Wc4A2tQ6WWY&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=7)
- [How to use modules in a pipeline](https://www.youtube.com/watch?v=tWvou0xj9wA&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=6)
- [Modules test data](https://www.youtube.com/watch?v=QXfAerydAT0&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=17)
- [Test modules](https://www.youtube.com/watch?v=pjhscKyWH74&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=18)
- [Where do I start writing my own DSL2 pipeline?!](https://youtu.be/Z_uPj7fAes8)
- [Development environments & workflows (Phil)](https://www.youtube.com/watch?v=XB96efweCLI&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=12)
- [Development environments & workflows (Maxime)](https://www.youtube.com/watch?v=OF55x-FT5WE&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=19)

# Hackathon groups

As the main objectives of this hackathon will be adding nf-core modules and rewriting existing pipelines to the new Nextflow language format, we have set up a few groups working on this task for different pipelines.

We will coordinate our work and the issues we are working on using a single GitHub [<i class="fab fa-github"></i> Project Board](https://github.com/orgs/nf-core/projects/20).

- **variant calling**
  - Working on adding modules for variant calling pipelines, mainly sarek, and raredisease.
  - Group lead: Maxime Garcia, Friederike Hanssen.
- **microbial genomics**
  - Working on adding modules for meta-omics pipelines, mainly ampliseq, mag, bacass.
  - Group lead: Daniel Straub.
- **epigenetics**
  - Working on adding modules for epigenetics pipelines, mainly chip-seq, methylseq, atacseq, HiC.
  - Group lead: Harshil Patel.
- **ancient DNA**
  - Working on adding modules for ancient DNA sequencing pipelines, mainly eager, coproid.
  - Group lead: James Fellows Yates.
- **single-cell**
  - Working on adding modules for single-cell pipelines, mainly scrnaseq, scflow.
  - Group lead: Gisela Gabernet.

# How we will work

We will be a lot of people working in parallel during this hackathon, so to stay organised we have a recommended workflow:

1. :speech_balloon: Most of the event will happen on [Gather town](https://www.gather.town/). You can chat there with your group to get an overview of what is going on.
2. <i class="fab fa-slack"></i> Join the `#hackathon-oct2021-public` Slack channel to stay up to date with the hackathon events.
3. <i class="fab fa-github"></i> Find a task to work on using the [GitHub Project Board](https://github.com/orgs/nf-core/projects/20).
    - If you have something you want to do that's not there, please make an issue (e.g. in the nf-core/modules repository if you are adding a new module) and add it to the board
4. :raising_hand: Assign yourself to the issue that you're currently working on (preferably one issue at a time)
    - This is so that multiple people don't accidentally work on the same task
5. :fast_forward: When you're done, make a pull-request with your changes. Link it to the issue so that the issue closes when merged.
6. :page_facing_up: Describe your work on the HackMD document for the project and tell the group! :tada:
    - [<i class="fas fa-file-alt"></i> Modules and Pipelines](https://hackmd.io/@nf-core/Bk42vnEHF/)
    - [<i class="fas fa-file-alt"></i> Tools](https://hackmd.io/@nf-core/S1vAu24BK/)
    - [<i class="fas fa-file-alt"></i> Website](https://hackmd.io/@nf-core/BkBHF2ErY)
7. :recycle: Repeat!

> The HackMD document is the easiest to forget, but please add something even if you think what you did was small -
> we will use it in the group check-outs for each day, to stay tuned between time-zones and also in the reporting after the event so it's important for us :bow:

# Schedule

We expect people to come and go during the hackathon due to diverse time zones.
Please just do whatever works best for you!

All locations described in the table below refer to places in the Hackathon gather.town
space.

The following schedule should display times in your local time zone:

<div class="table-responsive">
    <table class="table table-hover table-sm table-bordered">
        <thead>
            <tr>
                <th>Time</th>
                <th>Weds. 27 Oct, 2021</th>
                <th>Thu. 28 Oct, 2021</th>
                <th>Fri. 29 Oct, 2021</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td data-timestamp="1635379200" data-timeformat="HH:mm z">02:00</td>
                <td rowspan=4></td>
                <td>Check-out/in (Americas/Asia-Pacific)<br>
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Edmund Miller/Zhaowei Yang/Bhargava Morampalli</li>
                  </ul></td>
                <td>Check-out/in (Americas/Asia-Pacific)<br>
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Edmund Miller/Zhaowei Yang/Bhargava Morampalli</li>
                  </ul></td>
            </tr>
            <tr>
                <td data-timestamp="1635381000" data-timeformat="HH:mm z">02:30</td>
                <td rowspan=2>Pair-programming (Asia-Pacific)<br>
                    <ul class="small list-unstyled">
                        <li>Location: Mentor specific rooms</li>
                    </ul>
                </td>
                <td rowspan="3">Hack! <...></td>
            </tr>
            <tr>
                <td data-timestamp="1635382800" data-timeformat="HH:mm z">03:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635384600" data-timeformat="HH:mm z">03:30</td>
                <td rowspan="1">Hack! <...></td>
            </tr>
            <tr>
                <td data-timestamp="1635321600" data-timeformat="HH:mm z">10:00</td>
                <td>
                  Welcome
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Harshil Patel/Matthias Hörtenhuber</li>
                  </ul>
                </td>
                <td>
                  Check-out/in (Asia-Pacific/EMEA)
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Zhaowei Yang/Bhargava Morampalli/Matthias Hörtenhuber</li>
                  </ul>
                </td>
                <td>
                  Check-out/in (Asia-Pacific/EMEA)
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Zhaowei Yang/Bhargava Morampalli/Matthias Hörtenhuber</li>
                  </ul>
                </td>
            </tr>
            <tr>
                <td data-timestamp="1635323400" data-timeformat="HH:mm z">10:30</td>
                <td rowspan="2">Pair Programming (EMEA)<br>
                   <ul class="small list-unstyled">
                    <li>Location: Mentor specific rooms</li>
                  </ul></td>
                <td rowspan="2">Hack!</td>
                <td rowspan="2">Hack!</td>
            </tr>
            <tr>
                <td data-timestamp="1635325200"  data-timeformat="HH:mm z">11:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635327000" data-timeformat="HH:mm z">11:30</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li><li>Host(s): Susanna Marquez</li></ul</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li><li>Host(s): Susanna Marquez</li></ul</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li><li>Host(s): Susanna Marquez</li></ul</td>
            </tr>
            <tr>
                <td data-timestamp="1635327900" data-timeformat="HH:mm z">11:45</td>
                <td rowspan="3">Hack!</td>
                <td rowspan="3">Hack!</td>
                <td rowspan="3">Hack!</td>
            </tr>
            <tr>
                <td data-timestamp="1635328800" data-timeformat="HH:mm z">12:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635330600" data-timeformat="HH:mm z">12:30</td>
            </tr>
            <tr>
                <td data-timestamp="1635332400" data-timeformat="HH:mm z">13:00</td>
                <td>
                  Group catch-up<br>
                   <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Matthias Hörtenhuber</li>
                  </ul>
                </td>
                <td>
                  Group catch-up<br>
                   <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Matthias Hörtenhuber</li>
                  </ul>
                </td>
                <td>
                  Group catch-up<br>
                   <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Matthias Hörtenhuber</li>
                  </ul>
                </td>
            </tr>
            <tr>
                <td data-timestamp="1635334200" data-timeformat="HH:mm z">13:30</td>
                <td rowspan="3">Hack!</td>
                <td rowspan="3">Hack!</td>
                <td rowspan="3">Hack!</td>
            </tr>
            <tr>
                <td data-timestamp="1635336000" data-timeformat="HH:mm z">14:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635337800" data-timeformat="HH:mm z">14:30</td>
            </tr>
            <tr>
                <td data-timestamp="1635339600" data-timeformat="HH:mm z">15:00</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li></ul</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li></ul</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li></ul</td>
            </tr>
            <tr>
                <td data-timestamp="1635340500" data-timeformat="HH:mm z">15:15</td>
                <td rowspan="4">Hack!</td>
                <td rowspan="4">Hack!</td>
                <td rowspan="4">Hack!</td>
            </tr>
           <tr>
                <td data-timestamp="1635341400" data-timeformat="HH:mm z">15:30</td>
            </tr>
            <tr>
                <td data-timestamp="1635343200" data-timeformat="HH:mm z">16:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635345000" data-timeformat="HH:mm z">16:30</td>
                            </tr>
            <tr>
                <td data-timestamp="1635346800" data-timeformat="HH:mm z">17:00</td>
                <td>
                  Check-out/in (EMEA/Americas)
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Matthias Hörtenhuber/Edmund Miller</li>
                  </ul>
                </td>
                <td>
                  Check-out/in (EMEA/Americas)
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Matthias Hörtenhuber/Edmund Miller</li>
                  </ul>
                </td>
                <td>
                  Wrap Up
                  <ul class="small list-unstyled">
                    <li>Location: Lecture Theatre</li>
                    <li>Host(s): Harshil Patel, Matthias Hörtenhuber </li>
                  </ul>
                </td>
            </tr>
            <tr>
                <td data-timestamp="1635348600" data-timeformat="HH:mm z">17:30</td>
                <td rowspan="2">Pair Programming (Americas)<br>
                   <ul class="small list-unstyled">
                    <li>Location: Mentor specific rooms</li>
                  </ul></td>
                <td rowspan="2">Social Event<br>
                   <ul class="small list-unstyled">
                    <li>Location: Lounge</li>
                    <li>Host: Maxime Garcia / Susanna Marquez</ul></td>
                    <td rowspan="6"></td>
            </tr>
            <tr>
                <td data-timestamp="1635350400" data-timeformat="HH:mm z">18:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635352200" data-timeformat="HH:mm z">18:30</td>
                <td rowspan="2">Hack!</td>
                <td rowspan="2">Hack!</td>
            </tr>
            <tr>
                <td data-timestamp="1635354000" data-timeformat="HH:mm z">19:00</td>
            </tr>
            <tr>
                <td data-timestamp="1635355800" data-timeformat="HH:mm z">19:30</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li><li>Host(s): Susanna Marquez</li></ul</td>
                <td>Break<br><ul class="small list-unstyled"><li>Location: Cafeteria</li><li>Host(s): Susanna Marquez</li></ul</td>
            </tr>
            <tr>
                <td data-timestamp="1635356700" data-timeformat="HH:mm z">19:45</td>
                <td>Hack! <...></td>
                <td>Hack! <...></td>
            </tr>
        </tbody>
    </table>
</div>

## Social Activities

During the hackathon, we will have a few light-hearted fun and games!

- Breaks will happen in the dedicated Cafeteria room, for informal chatting and getting to know each other.
- Throughout the three days, we will once again be running a nf-core hackathon **bingo**! To join the game, you can go the following [link](https://nfcore-bingo.web.app/?game=nf-core-hackathon). Check the instructions at the bottom of the page.

    > <i class="fas fa-hand-paper"></i> Bingo! <https://nfcore-bingo.web.app/?game=nf-core-hackathon>

- In addition, we will be running an easter-egg (sock? 😉) hunt! There are 11 socks distributed around the gather.town world. Take screenshots of as many as you can find!
- Finally, during Thursday's social event (see schedule above), we will be running a short quiz!

All social activities are of course optional, but hope to see as many people joining in as possible :tada:
