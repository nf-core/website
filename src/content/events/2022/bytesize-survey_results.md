---
title: 'Bytesize: The Nextflow and nf-core community survey'
subtitle: Phil Ewels - Seqera Labs
type: talk
startDate: '2022-05-24'
startTime: '13:00+02:00'
endDate: '2022-05-24'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/DYdPNPVa4wc
locations:
  - name: Online
    links:
      - https://youtu.be/DYdPNPVa4wc
      - https://doi.org/10.6084/m9.figshare.19923410.v1
---

At the start of 2022, nf-core joined the Nextflow community to do a joint survey entitled _"The State of the Workflow 2022: Nextflow and nf-core community survey"_ (see the [announcement blog post](https://seqera.io/blog/the-state-of-the-workflow-the-2022-nextflow-and-nf-core-community-survey/)). This follows on from the [2021 edition](https://seqera.io/blog/state-of-nextflow-2021-results/) but with some additional questions, including some concerning nf-core community strength and engagement.

This week, Phil Ewels ([@ewels](https://github.com/ewels/)) will give an overview over the key results of the survey. If that's not enough excitement for you, Evan Floden ([@evanfloden](https://github.com/evanfloden) will wrap up with a live-draw to select the winners of the survey draw: one Apple Macbook Air M1 up for grabs alongiside ten sustainability-certified Seqera Labs merchandise packs.. 🤩 💚 💻 👕

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/DYdPNPVa4wct=1)
Thank you everybody for joining today's nf-core bytesize talk. Usually these short talks are about specific topics like how to develop code within the nf-core framework or about specific pipelines things like that. Today's a little bit different because we're going to go over the results of the Nextflow and nf-core community survey.

[0:20](https://youtu.be/DYdPNPVa4wct=20)
Those of you who are on Twitter and who are on the nf-core Slack will have seen me posting, requesting everybody to take a few minutes to fill in the survey back at the start of the year, around January / February time. Basically, it's something we're trying to do annually to really take a snapshot of all the different people using Nextflow, who everyone is, why everyone is using Nextflow, what works, what doesn't work, to try and prioritize development, and also really get a feel for what needs the most attention, both in the community and in the software.

[1:00](https://youtu.be/DYdPNPVa4wct=60)
So, apologies for that spam back at the start of the year, if you had that in multiple channels, but many thanks to everybody who filled it in. Those of you who followed through on clicking that link will have ended up on the Seqera labs webpage, which looked like this and you went through and followed the multi-step survey. One of the main reasons we want to do this survey is because Nextflow and nf-core are co-funded by a Chan Zuckerberg Initiative grant. This particular grant is called a diversity and inclusion grant through the "essential open source for science" program. And the focus of this particular grant that we're on is about trying to what the grant cycle is: improve the diversity and inclusion both geographically and through every other metric.

[1:53](https://youtu.be/DYdPNPVa4wct=113)
In order to track whether we're doing a good job, we need some metrics. It's very difficult to track this, but one of the things we want to use is this survey, basically. And so by doing it pretty early on within the scope of this two year grant we're hoping we can track improvement over the next two years, time will tell. It's really important for us within the context of community growth and funding. Let's dig in.

[2:21](https://youtu.be/DYdPNPVa4wct=141)
Those of you who are active on Twitter may have noticed that a Seqera labs tweet went out a couple of hours ago, there is a blog post on Seqera labs website, all about this with the infographics. You can find all of this information dig into it yourself. If you haven't done already. I'm going to go through some of the key conclusions in this talk. And I'm also going to put out a few additional statistics, which didn't make it into the infographic, just so that you don't feel like I'm repeating myself completely.

[2:49](https://youtu.be/DYdPNPVa4wct=169)
Let's start off with some community demographics. Firstly, as hopefully we already knew and hoped, we are very global community, which I love. The majority of users are based in the US and in the UK and certainly in Europe. That's fairly inevitable from our origins of the community and also mirrors the density of people working in the field in bioinformatics. But there was I think 36 different countries in the respondents list, which is fantastic. I'm sure that's up a lot since the last few years. That's really nice to see. We're increasingly spreading around the world. And let's see if we can push these numbers up and make that map go even more blue then next year.

[3:37](https://youtu.be/DYdPNPVa4wct=217)
The majority of people, we asked what your favourite primary language for reading and writing was. Most of you picked English, which is not a surprise. But there's quite a lot of people speaking other languages as well. We have a pretty terrible gender equality. I'm not sure it's just our fault. I think it's probably indicative of a wider issue, but that's definitely something that could, of course, be improved if anyone has any ideas. And yeah, there's a pretty wide range of people, lots of early stage researchers using Nextflow and lots of people well into their career as well. It's really nice to see these kinds of things and get a feel for who everyone is.

[4:19](https://youtu.be/DYdPNPVa4wct=259)
So I said there's lots of languages. You can see English is up there at the top, but there's a long, long tail. And interestingly, a lot of people bundled into that "other" category there. So, again, there's a nod to how diverse our community is already. And this is really useful, for example, if we want to prioritise any efforts to translate material. We know which language is the most important to our community.

[4:45](https://youtu.be/DYdPNPVa4wct=285)
Digging into a bit more of what it is that everybody does, very similar to last year, the majority of people who filled in the survey are classed themselves as bioinformaticians. Few other people with different categories, job categories is always difficult. You can look into the others category of people, some identity crisis issues going on there. But most people are bioinformaticians working with biological datasets. Lots of people within academia and research, but also a lot of people in biotech startups, especially that seems to be growing since last year, and pharma and clinical work. That's really interesting to see as Nextflow matures and gets more heavily adopted, it's branching out of academia a little bit into the wider community.

[5:35](https://youtu.be/DYdPNPVa4wct=335)
Lots of people who filled in the survey have only recently started using Nextflow, which is really interesting, still under a year for the majority of people who filled in the survey. Welcome, all of you. Even though some of the statistics came out similar to last year, we're actually looking at a lot of people here who are new. And I think that's fantastic. It shows we're still, I haven't saturated the market by any means. There's still lots of people who don't know about Nextflow and lots of people joining the community all the time.

[6:08](https://youtu.be/DYdPNPVa4wct=)
Generally, you're a very happy bunch. Everybody likes Nextflow, which is good. Maybe there's a bias in who fills in the survey here, but generally everyone seems to say that they're very happy with Nextflow and with the community. The vast majority of you would recommend - and I believe do recommend - Nextflow to your colleagues. And that's actually slightly better than last year. An even slightly better satisfaction rate, which is never a bad thing. Always room to improve at the top.

[6:40](https://youtu.be/DYdPNPVa4wct=400)
Something that didn't make it into a blog post, but I think is one of the more interesting parts of the survey, is those of you who felt frustrated with Nextflow. It's not a complete paradox here. I think it's fine to say you're satisfied with Nextflow, but you are occasionally feeling frustrated with it. That's natural with any programming language. And so if you've ever felt like this, you're not alone. Most of us have at times felt frustrated with Nextflow.

[7:05](https://youtu.be/DYdPNPVa4wct=425)
I dug into that a little bit and started reading 300-and-something free text responses here about why all of you have felt frustrated. The common themes that jumped out to me were familiar to many of you, I'm sure, the fact that Nextflow works with Groovy, which is not one of the mainstream languages for bioinformatics people. And a lot of people say they often struggle to interpret what the error messages mean. It's people saying it can be quite difficult to get into Nextflow and nf-core as a steep learning curve. And a few people saying that the community is so active, things are moving so fast, it was difficult to stay up to date, which is a double edged sword there. There's lots of activity, which is great, but it can be difficult to keep up. Just for those of you who filled in this question, know that we hear you. These are all things we're aware of within Nextflow and nf-core and things that we're always trying to improve on.

[7:59](https://youtu.be/DYdPNPVa4wct=479)
I went a bit further, just because I could, and threw together a quick word cloud here for all the things that annoy you, just as a form of venting, I guess. But it can't be all that bad because you're all really happy. A couple more questions here. People asking why you're running Nextflow. And the first two categories, people running and writing their own workflows, basically doing analysis for themselves. The next two categories are people running and building workflows for other people. Bioinformatics core groups and things like that. And a handful of you building larger systems that include Nextflow.

[8:43](https://youtu.be/DYdPNPVa4wct=523)
When it comes to the workflows you're using, lots of people building their own workflows, which is of course expected, but a fantastic number of you using nf-core workflows. Now, again, in fairness, there's probably some bias here. We've pushed out this survey through nf-core channels amongst others. I would sort of hope that at least some of you were using nf-core pipelines. But still, it's fantastic to see so many people responding that they are using nf-core workflows on a regular basis. And this is a really valuable resource.

[9:13](https://youtu.be/DYdPNPVa4wct=553)
But you're also quite promiscuous. It's not just Nextflow. Over half of users are using more than one workflow tool, which I was quite surprised by this result. Lots of you are using Snakemake and Nextflow, Galaxy, CWL, whatever you need to get the job done. Just because you're very happy with Nextflow doesn't mean you're blind to all the alternatives. And that's no bad thing. It's good to have some competition and cross fertilization of ideas.

[9:42](https://youtu.be/DYdPNPVa4wct=582)
When it comes to where you're running compute environments that everyone uses, the majority, just like last year, are still running on HPCs, on-premises clusters and also single computers. And that's not something we expect to change massively in the near future. But there is an uptick in the number of people using cloud. If you compare to last year, the categories are pretty much the same, but basically there's a bit of an increase in the people who are using private clouds, especially. Anyway, decreasing the number of people using HPCs.

[10:22](https://youtu.be/DYdPNPVa4wct=622)
For those of you who are running with HPCs, the majority use SLURM. That's definitely the most common scheduler, followed by Grid Engine, and that's again similar to last year. And we see that in the community on Slack, people posting questions, lots of people using SLURM. Quite a lot of people using public cloud today, and quite a lot of people planning to move towards the cloud, lots of people looking in that direction. And when we compare the different types of cloud to different public clouds available, AWS, Amazon is by far the most popular, but Azure has a climbing rank there. And again, if you break this down by where people are answering from, lots of people in academia working with clusters and public cloud is super popular within the private sector, which maybe is not that surprising. But yeah, up to 77% of people within private sector are using the public cloud now.

[11:23](https://youtu.be/DYdPNPVa4wct=683)
Last year we asked about Kubernetes, a bit of a hot topic for those of you who know about it. Lots of people, a small number of people, I think it was about 8% last year, who are already using Kubernetes. Lots of people saying you're planning to use Kubernetes in the future. We're curious to see if anything changes in the year. It hasn't. There's pretty much exactly the same number of people saying that they are actively using Kubernetes today across the various different Kubernetes solutions. But again, lots of people saying that they're interested in it. Moving forward, we'll see if anything changes there in the future.

[12:07](https://youtu.be/DYdPNPVa4wct=727)
We asked a bit about the different traits, the reasons that you're use Nextflow, what things do you find important when you're choosing which tool to use. And the winner of this category was definitely documentation. I'm totally with you on this one. I've got a soft spot for trying to put together documentation for tools. Documentation and performance are the two outstanding categories here. And when we asked you, OK, within documentation and learning materials, what do you use the most, what's most useful? The official documentation, the reference documentation came top and nf-core came a close second. That's great to see. Everyone is using the documentation that we've been building. And this is definitely a hot topic for us right now. We've got lots of room for improvement with documentation. That's good to see.

[13:08](https://youtu.be/DYdPNPVa4wct=788)
This survey went out just before the Nextflow slack went live. Everyone was still using the Nextflow gitter at that point. But we asked about the nf-core slack. And the vast majority of people who responded knew about and used the nf-core slack, which is great. Not only is the nf-core community building pipelines and standards, but it's also a big support channel. And lots of people are also feeding back into the community. Nearly 30 percent of people are contributing back to these nf-core community pipelines. It's really good to see. It's a two way street there.

[13:46](https://youtu.be/DYdPNPVa4wct=826)
Apparently no one really wants a graphical interface for their workflow manager. It's fair enough. But you can see documentation really stands out as being very important to lots of people. And then also we're thinking about integrations and this is more tooling, it's quite detailed here. But a lot of you want to be able to optimize computational resources, which makes sense. Unit tests was a popular category here. This is great. This is real fodder here for feeding into the Nextflow development process, and nf-core, to really prioritize which topics need to be tackled.

[14:33](https://youtu.be/DYdPNPVa4wct=873)
Right. I'm going to wrap up there. You can go and look into this in more detail yourself and make your own conclusions on the Seqera website. And we've got the blog post which went live this morning so you can click through to that. If you have any questions, just quickly check the questions now. Before we go on. No? OK, in that case, I will hand over at this point to Evan, CEO of Seqera, and Evan is going to share with us a live draw for the present surprises because there was definitely an ulterior motive for some of you to fill in this survey. And with that, I'll stop sharing and I'll pass over to Evan.

[15:21](https://youtu.be/DYdPNPVa4wct=921)
(Evan) Thanks a lot, Phil. This should be pretty short here. But what we have is part of the survey we had. Around 10 people, 10 prizes, which were set up for our Nextflow and Seqera Labs merchandise packs. As part of that we took all of the names, and we placed them into this big sort of circular prize drawer here. And from there we drew 10 people. We just did this earlier on this morning, mostly because it takes some time, but if you want to watch it you can go through the whole thing. The winners of the merchandise packs that had as part of that are Jacob, Stefano, Susanna, Anca, Li Z, Yuk K, Niclas, Adam, Avinash, and Chela. We'll reach out to all of you, send you out a link where you can get this. It's got a hoodie in there. We've got t-shirts, hats, and some cool Nextflow socks as well. Thanks to everyone for doing that.

The next part of it is around the prize for the Mac. As part of that we're going to do that live and I'm going to share my screen here and hopefully be able to do this. Again this is the similar thing where everyone's names be entered into this, and we will draw the winner. Hopefully you can share a part of my screen here to do that. Okay wasn't too far off. So, this will take about 10 seconds or so to run through, and we will have the winner. Let's start that off now. Let's spin around. And the winner of the Mac M1 for the prize is Michael H from USA. Congratulations Michael we'll reach out to you and send you through an email with the information and go from there. That's it. I think thanks so much to Phil for bringing that presentation together.

As I say, if you want to go reach out and have a look at the blog. We are going to try and do these more regularly. This is the second year of this and a lot of this information becomes more useful the more as it's worth. I'm not sure if for some reason my video is off but here I am. I'm a real person. Thanks everyone for joining. Say, read the blog, reach out to us if you've got any questions, always appreciate everyone's work and yeah thanks so much.

[17:45](https://youtu.be/DYdPNPVa4wct=1045)
(host) Thank you very much. Are there any questions for Phil or Evan from the audience? I don't think so. Then I would like to thank both of you, and also the Chan Zuckerberg Initiative for funding, of course, and I'm going to stop the recording now, but be aware that you can always ask more about today's topic in the bytesize channel on Slack. And yeah, contact us if you have any other questions.

</details>
