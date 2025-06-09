---
title: "Bytesize: nf-core hackathons"
subtitle: Friederike Hanssen, QBiC
type: talk
startDate: "2024-01-30"
startTime: "13:00+01:00"
endDate: "2024-01-30"
endTime: "13:30+01:00"
youtubeEmbed: https://youtu.be/nmp0pXL7GqQ
locations:
  - name: Online
    links:
      - https://youtu.be/nmp0pXL7GqQ
---

Curious what we do at nf-core hackathons? Ever wondered if you should sign up? In this weeks bytesize, Friederike Hansen ([@FriederikeHanssen](https://github.com/FriederikeHanssen)) will have all the answers for these questions and more!

<details markdown="1"><summary>Video transcription</summary>

:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/nmp0pXL7GqQ&t=1)
Welcome, everyone, to this week's bytesize talk. With us today is Rike Hansen from the QBiC in Tübingen, and she's going to talk about everything you need to know about hackathons. Keeping in mind, there's one coming up very soon. Off to you, Rike.

[0:16](https://youtu.be/nmp0pXL7GqQ&t=16)
Hey, thank you. Or at least I'll try to talk about everything relating hackathons. To get started, maybe quickly, what do we mean by hackathon? I don't know what the official definition of hackathon is, but what we mean is we just want to meet up for a few days with other like-minded people and spend some time coding together, maybe actually close some issues that have been hanging out for a long time and just enjoy each other's company and work together. Who is this aimed at hugely? We aim hackathons at people that are experienced in Nextflow and are interested in contributing to nf-core pipelines, to various nf-core infrastructure things, or also to Nextflow plugins. If you're completely new to Nextflow, the hackathon might not be quite the right fit for you yet, but instead we have run several community trainings a year, usually also around the hackathon before, that you can attend. All the material for this is also online and the recordings of it are on YouTube. You can go back, rewatch it, or also do it whenever you have time to do so.

[1:30](https://youtu.be/nmp0pXL7GqQ&t=90)
We have two setups for hackathons. We do some in person where most people meet and then we do some online. We do the in-person one, usually in Barcelona in fall, roughly. Then we also recently started having one in Boston. The other one that we do is online. For this, we use Gathertown, just like a little platform where you can log in, you can have your little avatar that you see there. You can walk around, you can easily interact with people. You can go to enclosed spaces to communicate with a bunch of people, but don't bother anybody else. It's really easy to use to just get this hackathon and group feeling while being online. It also has a bunch of fun little things like dogs you can pet and Go-cards you can race around. We've been using gather for quite some time and it's worked really well. There are a couple of things to keep in mind. You need to use the email address that you registered with for the hackathon. If Gathertown is sounds a little bit complicated or you don't really know how to use it. James made an entire bytesize just about gather that you can check out.

[2:49](https://youtu.be/nmp0pXL7GqQ&t=169)
Online hackathons are also usually once a year and they tend to be in March, let's say spring. Last year we extended the online hackathons to having distributed local sites. We encourage people at the institute, at their workplace to have small gatherings of people that are around to mix a little bit this online only, but hackathon together with the in-person feeling of it. These are pictures from last year and by making the slide, this was probably the slide I had most fun with going back to and finding all the group pictures from the hackathon in March last year. If you want to host one of these distributed sites, it's very low stakes and low effort to make it as easy as possible. Essentially the only two things you need to do if you want to host one is you need to book a room and you need to add yourself to the website. You can see here already a bunch of people have registered for the upcoming hackathon. Then when people sign up for the hackathon, they can find your site and register there up until the maximum number of people has reached. Then you just program together. We do ask people to be online on gather while they meet, to stay in touch with everyone else, since the majority of people will not be in the room with you. So, don't forget to bring your headphones to not bother anybody.

[4:10](https://youtu.be/nmp0pXL7GqQ&t=250)
How do we normally do hackathons? We subdivide ourselves and topics, formerly groups. Typically we have a topic or a group that centers around pipelines. People that want to work on existing or new nf-core pipelines, then we have a topic that centers around modules and sub workflows, topic on infrastructure, nf-core tools and website development. We tend to have a few more depending on what we're doing at the time. Recently we had topics focused around nf-tests or Nextflow plugins or documentation. You are not at all bound by it, it just helps us organize a little bit. Throughout the hackathon, you can switch around as you want. You might develop modules to be used in your pipeline that you then like to add on. We just use these topics to organize, but you can do whatever you want and just move around as you please.

[5:27](https://youtu.be/nmp0pXL7GqQ&t=327)
We communicate as always in nf-core via Slack. We set up specific hackathon channels. Here I already made up some that we will probably have for the upcoming hackathons. We have a central channel where general stuff like announcements or wrap-up is happening. These things are posted. Then for each topic, we have a channel to find people that are working on the same thing or to get help or just talk to people there. We organize all of our issues and features in GitHub project boards. For each hackathon, we create a new one. If you go to the nf-core GitHub site, you will find this tag project. Then there you will find the link to the hackathon project board. There we collect any issues or features that people want to work on. You can also add your issue if you want to work on this during the hackathon. It helps us to organize, it helps other people find issues if they're looking for work. It helps to make sure that not two people work on the same thing. It's really nice to have. Also it helps us track what needs reviews, what has been done, and just the current state of things. If you can't find it, we also have a website for every hackathon where the Slack group is linked, all the topics are listed, and also the link to the project board is listed. Here at the bottom, you can see it in gray.

[7:01](https://youtu.be/nmp0pXL7GqQ&t=421)
During the hackathon, we have all of these issue cards, these are the same issues that you have in your pipelines or modules, they're just linked there and we can use them there. If you pick an issue to work on on the project board, don't forget to assign yourself just so everybody else knows that somebody is working on it, or you also can more easily find people to help you out on these issues and then also update the status. You see the little drop down menu that says "to do", "in progress", "done", and "in review". Just by updating the status, it's easier to keep track.

[7:39](https://youtu.be/nmp0pXL7GqQ&t=459)
One of the probably the most important parts actually is reviews. These are also the ones that usually end up being quite a bottleneck. We need to review everything, all the pull requests before merging them in. We need people to do this. Review as you go. Also in some groups, we have tried out review buddies. For example, two pipeline groups, RNA-Seq and raredisease, so that they exchange reviews with each other to get them done more quickly. Also don't forget to add the ready for review in your status and drop the link on the request review channel on Slack just to raise awareness and also return reviews. If you just worked on a module and open a PR, maybe take five minutes to pick another PR to review. If you're unsure if you should give an approving review, you can also just comment on it. This is also really helpful for other people.

[8:42](https://youtu.be/nmp0pXL7GqQ&t=522)
Throughout the Hackathon day, you will be bothered by your topic leaders to fill out some progress slides. At the end of each Hackathon day, we have some up slides that just give a very high level overview of the day, what happened, which modules were worked on here, for example, what that Hackathon, we were particularly talking about nf-tests. We highlighted this or other interesting things like a new sub workflow was added. It is there to update everyone else and also to celebrate ourselves a little bit for the achievements.

[9:39](https://youtu.be/nmp0pXL7GqQ&t=579)
To sum up how to contribute at a Hackathon: Go on Slack chat with your group, then find an issue and assign yourself and start working on it. Open a PR and request a review. Maybe while you're waiting for the review, review one back and then update the progress slides. Then don't forget to celebrate, get up, get a cup of coffee, get some snacks, take a walk and then start over again. All the information for all the Hackathons is also always linked on the website. For each Hackathon, you will find an individual site. Here's the March one already linked.

[10:02](https://youtu.be/nmp0pXL7GqQ&t=602)
A couple of things before I finish, if you're very familiar with Nextflow, but you are not so familiar with nf-core yet, we have a bunch of bytesize talks that cover individual aspects, like how do we do testing, linting, how do people have their code environments set up, what is nf-test, all these things that could be very useful. Maybe check out those. If you're planning to work on a new pipeline, be sure to propose it in the new pipelines channel before the Hackathon. Whenever somebody has a new pipeline, they get proposed and then we talk about it. Does it fit into another pipeline? Can it be a standalone pipeline? What should the name be? And it takes a little bit of time and if you propose a new pipeline at a Hackathon, you might not actually be able to work on it because people didn't have time to review it or discussions took too long. If you want to work on a new pipeline at a Hackathon, come to the new pipelines channel and write about it there. Also maybe when you sign up for the Hackathon, consider a little bit what you want to work on. It's nice to get in touch with others beforehand and find some collaborators on something and then you can start right away on Monday morning. Last but not least, I hope this goes without saying, but we have a code of conduct that says, don't be mean, be inclusive, be friendly, and everybody who attends one of our events needs to adhere to it.

[11:30](https://youtu.be/nmp0pXL7GqQ&t=690)
As Fran mentioned, we have a Hackathon coming up in March from the 18th to the 20th and we are still looking for a lot of local site organizers. If you're thinking about it, I can definitely recommend it. We did it last year and it wasn't a lot of effort and we had a lot of fun to do it. If you have questions about it, I'm also happy to answer anything there. Okay. Thank you very much.

[11:54](https://youtu.be/nmp0pXL7GqQ&t=714)
(host) Thank you so much. I have now allowed everyone to unmute themselves to start their video. If there are any questions, please ask away.

(speaker) Any in the chat?

(host) There's a lot of thank yous in the chat, but I have a question. I actually have a question.

(question) If you are interested in working with Nextflow and nf-core, but you don't actually have your own pipeline, you just learned Nextflow with one of the courses, what is it you can do at a Hackathon?

(speaker) So one thing that I think always works really well is pick some of the modules and contribute those. They're very self-contained packages usually. You have one module, it's very achievable to add this one tool, I think, minus a couple of difficult points like test data or so. For this, you don't need any pipeline. The modules are completely independent of pipelines and you can always contribute there. But also if you want to contribute to a pipeline, a Hackathon is a great way to meet the developers, to get in touch with them and work together with them as well.

[13:14](https://youtu.be/nmp0pXL7GqQ&t=794)
(question) Thank you very much. Yeah. I have a question. If we have like an idea to work on something, what's the proposed procedure to actually start on to work on like a certain component of a pipeline that we feel is missing or something like that?

(speaker) Like maybe adding a new feature to an existing-

(question cont) Yes, for example, yes.

(answer) Add an issue to the GitHub repository of the pipeline, maybe talk to the people who developed the pipeline if something is already going on and then it will be added to the issue board and then you can just work on it. If it for some reason didn't end up on the project board for the Hackathon, don't worry. You can still work on it and we'll add it later.

(question cont) Thank you very much

[14:02](https://youtu.be/nmp0pXL7GqQ&t=842)
(host) There is a question in the chat. Can we work on nf-core style workflows, for example, nf-core template-based workflows using nf-core modules?

(answer) Sure, nobody's keeping you from anything. I guess we can say that any nf-core related work is very welcome at the Hackathons. It often doesn't have to fit in any of the specific topics that we choose. You can work on anything. It's more that we try to organize it a little bit to help people, but if you have something you want to work on, go for it. If it isn't on the project board, don't worry about it. Just add it later. Or if it doesn't fit there, also good. Maybe you found something completely new to work on that you just came up with back then like a completely new topic. It's all good. All very welcome.

[14:58](https://youtu.be/nmp0pXL7GqQ&t=898)
(question) Hello. I'm fairly beginner. I've played around with nf-core. I've written my own very simple, like read QC pipelines in Nextflow just to see how things are working. Do you think a Hackathon would be suitable for someone like me?

(answer) I think so. I mean, if you already know Nextflow a little bit and you know nf-core a little bit, I think it's perfectly fine. We just want to make clear that we don't have any training sessions at the Hackathon and we don't have typically people there that are dedicated to train other people. But if you have a few components you want to work on or just chat with people, then that's good, I think.

(question cont.) Cool. Just to follow up, do I need to be added to the nf-core GitHub to be able to?

(answer cont) Yes. Very good point. Typically, before you start a Hackathon, I meant to mention this, before you join the Hackathon, make sure you sign up to Slack and to the GitHub repository. Typically, there's also, an email or maybe it's in the signup sheet that you need to do it. But yes, there's two very important points. Thanks.

(host) I also want to clarify, even though we're not there to have the training events or anything if you're a beginner and you have questions, of course, ask those questions straight away. Don't hesitate. We're happy to answer any questions anyway, I mean, even outside of Hackathons. But there you might actually find people in the room that can immediately answer the questions.

(question cont.) That's actually awesome. Excellent.

[16:45](https://youtu.be/nmp0pXL7GqQ&t=965)
(question) Can infrastructure and DevOps professionals participate in the Hackathon or is this Hackathon focused solely on developing new code?

(answer) No, they do and they have participated on it or developing, I don't know, and maybe you mean like GitHub actions or something for deployment and so on. We have a bunch of people who work on this. I think actually last Hackathon, we had a group that was dedicated to DevOps.

[17:13](https://youtu.be/nmp0pXL7GqQ&t=1033)
(host) Are there any more questions? No, I think we covered it all. Then I would like to thank you, Rike, for this great presentation and everyone in the audience for listening and as usual, the Chan Zuckerberg Initiative for funding our bytesize talks. If you're interested, don't forget to sign up for the Hackathon or maybe even to host a site. Thank you.

</details>
