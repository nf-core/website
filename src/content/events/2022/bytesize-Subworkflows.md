---
title: 'Bytesize: Subworkflows'
subtitle: Maxime Garcia, Seqera labs
type: talk
start_date: '2022-11-22'
start_time: '13:00+01:00'
end_date: '2022-11-22'
end_time: '13:30+01:00'
youtube_embed: https://www.youtube.com/watch?v=-vHAXsuYQhE
location_url:
  - https://doi.org/10.6084/m9.figshare.21617526.v1
  - https://www.youtube.com/watch?v=-vHAXsuYQhE
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: Subworkflows

This week, Maxime Garcia ([@maxulysse](https://github.com/maxulysse)) will share with us all there is to know (well, some of it) about the brand-new nf-core subworkflows!

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=1)
(host) Hello, everyone. My name is Franziska Bonath. I'm today's host of the bytesize talk, and with me is Maxime, newly member of the Seqera team. He is here today to talk about nf-core subworkflows. Very interesting.

[0:18](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=18)
Indeed. Hello, everyone. Let me share my screen. I'm sharing my full screen because I'm trying some demo as well at some point, so I think it's a bad idea, but I need to try that. Full-size screen. Hello, everyone. My name is Maxime Garcia. I'm working at Sequera Labs in the Scientific Development team. I'm still working in Stockholm, but remotely from Barcelona, but yeah, that's Seqera, so it's fine. Come join us.

[0:58](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=58)
I'm going to talk about subworkflows. Basically, what is new with the subworkflows in nf-core? What is our plan, currently and maybe in the short, long, medium term, and some demo time. Tiny disclaimer to start, I think it's always important. More or less, what I say are more or less like my own takes of what the community is doing. Other developers might follow my ideas. Other developers might have other ideas. But I think it's good because that's how we are forging the best practice. It might and it probably will evolve some of the logic, some of the syntax, some of the stuff. But I think it's good that we try stuff and we figure out what is the best way to do stuff.

[1:48](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=108)
All of the new stuff. I think the most important part is that we have now subworkflows in a specific folder in the nf-core module repository. We can have a look there. If we look at the repo, we have a subworkflows folder. And in the subworkflows folder, we have an nf-core folder. And in that folder, we have all of the current subworkflows for nf-core.

[2:20](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=140)
I think the most important thing that we did during - not the most important thing, but for me, the stuff that had the most impact on my work on the pipeline - was actually the naming convention that we had, which is basically explained all there. We want to have a naming convention because that way it's much easier to understand directly what a subworkflows is doing or not.

[2:55](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=175)
Otherwise we have a lot of documentation here on this page, which is... I'm at the bottom of the page... Which is DSL2 subworkflows that's in the doc contributing... DSL2 subworkflows. A lot of the logic is inherited from the module. We still keep the same terminology, which I think is super important. Remember in nf-core a module will be just the atomic process and the subworkflow is a chain of modules. All of the logic is pretty similar to all that. All of the underlying logic.

[3:42](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=222)
We also have a lot of documentation for that. What we did new is some new tools commands. We don't have everything there yet, so I will... Let's finish the presentation before we actually start the demo. A new command is installed, so I will show that in the demo time as well.

[4:05](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=242)
The plan. We do have a plan of course. At least I have a plan myself, which is translating all of the local subworkflows that we had in Sarek and putting all of that into nf-core because that way I think it's very good for the community and I am hoping also to help other pipelines by doing that, and to convert more and more local subworkflows into nf-core. That way I'm pretty sure we can find a proper logic to be smarter and to do smart things. It's a bit redundant, but that's my plan. And I think with that we could find some new ways to do stuff for a new pipeline.

[4:50](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=290)
One example is that in Sarek we use Freebase. It's one of the tools that we use. We use many tools, but at least I know that Freebase is used by other pipelines. And I'm pretty sure that we're doing some stuff in Sarek. We have a spread and gather solution that can speed things up. we are trying to import that into the nf-core module, and I'm pretty sure once we have that in the nf-core module, then other pipelines might be able to import this subworkflow, and I'm pretty sure that will be a huge gain for the whole community. I'm really looking forward.

[5:26](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=326)
And now let's go into demo time. I'm going to demo how to install a subworkflow. I'm going to install a subworkflow that I just created yesterday. Actually, I just ported a local subworkflow from Sarek into a nf-core module. I'm just going to do that. First I already did that, I installed the dev version of tools. This part is the most important, I think, dev. I did that already. Now I'm going to my local repo. This is my own fork. I'm just creating a new branch. That's very simple. Now we have the new nf-core tool. There's this command `subworkflows` that the infrastructure team, so Matthias, Julia, and everyone else, I don't know who else is involved in the infrastructure team, but they did a pretty good job with all that. And I always like what they do. It looks so fancy, what they're doing.

[6:43](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=403)
No command subworkflow? Oh, yes, maybe just without the "S", subworkflow, no? subworkflow. And yes, the S was there. Don't misspell stuff! For the pipeline, we have `info`, `install`, `list`, `remove`, and `update` to develop new subworkflows that will be very similar to the same command that we have for the modules. I will not show that, but I will show all of that. Let's try `nf-core subworkflow info`. I want to have an info about the new subworkflow you want to install. Is the subworkflow locally installed? No, because I want to install it. Please select a subworkflow. I want to select the "vcf_annotate_snpeff", and I have some nice information about all that. It does perform annotation, I mean, snpEff, and then bgzip plus tabix index resulting VCF file. That's perfect. We do need a metamap vcf, a version of the snpEff database, an optional path to the root cache folder for snpEff, and then we have output. Compressed vcf file plus tabix index, html report, and of course the version.

[8:13](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=433)
What were the other commands that we could see? `Install`, `list`, `remove`, `update`. Let's check `list`. I want to list. List local. No local, that's some logic. And `list remote`, that's the same one. Then let's go for `install`. Wait, before I actually install, let's remove the one that we had. Git remove subworkflow local vcf_annotate_snpeff. I removed my... I removed my local version of this subworkflow. I will now install the new version. that's `subworkflows install`. And I want "vcf annotated". I could copy from there, but I want to try out what's happening if I don't do anything. This is so fancy! I love that! I am a big fan of auto-completion and stuff that is doing that. That was super fast, so it's done.

[9:42](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=582)
Use the following statement to include the subworkflow. I will grab that copy. I launch my code. I have a subworkflow, I have that local. This is my meta subworkflow where I do everything. I'm just copy pasting this new command here. Let's just align everything well. That looks good... Annotate all... Here we'll use this new subworkflow which is located there. That's nice. I don't have my "vcf_annotate_snpeff" anymore because it's now there. That's wonderful. That's what we expect.

[10:36](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=636)
Where are you with the kit? So we deleted some file, we added some new file, and we modified some file. Let's add the new file. We move to file. The meta.jml is different enough so it's like a new... It doesn't register as the renaming. The script itself for the subworkflow is exactly the same which makes sense because I created it yesterday and I basically copy pasted everything. What is happening in the module.json file? Get this module. This is new. Installed by module, installed by module. That's interesting. It's just looking, I like to check everything. I think it's important. So let's add this new file. Let's commit everything. Lets push. Let's create the pull request. We want to do that in dev. Let's replace subworkflow. Create the pull request. That looks good.

[12:41](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=761)
I'm thinking there is just one lesson that I need to do but this is very specific to Sarek. Yes, I need to change the path to the file here. We are doing pytest with tags and we are watching if some of the files are being changed or not from one PR to another. And then we are triggering the test just on that. For that, because the path is not the same anymore, I just update the path. This is done. Let's commit that as well. Let's push. I'm hoping that we are done with this pull request. Yes, we can see that it was failing before. I'm pretty sure because of the tests that were failing. Now everything is triggered. I can see pytest workflow is being triggered at the moment. I'm guessing once it triggered, it will figure out which test it has to run or not, but that's something else. I think that's good for that.

[14:20](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=860)
Let's go back to here. I paste the rest of the history here in my slides, and I will share my slides after this talk. I think now it's time to thank everyone and to go for the questions. This are the institutes that are participating in nf-core. I really need to update that slide because I think it's already one year old and I'm pretty sure we have like more people now. Same with the contributor, but I really want to thank everyone that is contributing to nf-core because it's a community and that's a community effort and without everyone else we wouldn't do anything. If you have any questions, please ask them because that was mainly just a demo and that was fairly simple. I'm pretty sure people have more questions.

[15:14](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=914)
(host) Thank you very much. You're now able to unmute yourself. If you have questions, either put them in the chat or ask them straight away. I think I saw some questions. It's not a question, it was a comment. Someone is being very happy that there are already 24 subworkflows.
(speaker) Yes, because we started the sub workflow at the hackathon properly. That was when? Last month or two months ago?
(host) Last month.
(speaker) Yeah, so 24 just in a month. That's good. I'm pretty sure we'll have more and more coming. And I know that Mathias is working on adding the command line help for nf-core tools soon. I'm guessing we're waiting for release of tools for that. John, do you have a question?

[16:14](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=974)
(question) Yes. Can you hear me? Thanks. Very interesting talk. I'm quite new to this, but I use Nextflow and I am also a little bit used to nf-core. But this thing about subworkflows, is this specific to nf-core or is it something that can apply to other Nextflow pipelines?
(answer) This is something that can be applied to any pipeline. We developed that first with nf-core in mind, like the module, but then every module, like in Nextflow, everything can be a module. Every process can be a module. Every chain of process can be a module. Even the workflow itself can be a module. You can import whatever you want, however you want. Definitely what we are creating here with nf-core, like this subworkflow stuff, it can be used in the broader Nextflow community without any issue.
(question continued) Okay. Thanks.

[17:15](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=1055)
(question) I also have a question or maybe comment also. That was a great presentation. However, I was going to ask, maybe my first comment is similar to what John just said. The presentation sounded more like subworkflows were nf-core things instead of a Nextflow thing. I think that's why he was asking that question about whether subworkflows were nf-core or Nextflow. My other question is that, what's the naming convention for subworkflows in nf-core? Is it like the first word is a verb followed by the names for the tools that you are chaining together in that subworkflow? Because I noticed some pattern like that, but maybe I'm wrong.
(answer) Yes, we have this convention. It's definitely like an nf-core thing only. I'm guessing like other people that develop stuff might want to follow the convention as well. I'm happy to talk more about that. But I think we have this convention. I think it's the input file type, which is the first. Then it should be a verb and then the list of the tools that are used. For example, like in that case, what we were doing with this subworkflow that I just added, it was vcf underscore annotate underscore snpeff.
(question continued) Yeah, thanks.

[18:51](https://www.youtube.com/watch?v=-vHAXsuYQhE&t=1131)
(host) Thank you very much. Are there any more questions? It doesn't seem so. If you have more questions, as usual, you can go to Slack, either in the bytesize channel or there's actually a channel also for subworkflows? (speaker) Yes, there is a channel for sub workflows. A channel for tools as well, obviously.
(host) Obviously. Or you can directly ask Maxime. Otherwise, I would like to thank Maxime for the talk and, as usual, for funding the Chan Zuckerberg Initiative. And you all for listening. Thank you very much.

</details>
