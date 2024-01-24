---
title: 'Bytesize 1: Introduction to nf-core'
subtitle: 'nf-core/bytesize: Bite-sized talks, (giga)byte-sized science!'
type: talk
startDate: '2021-02-02'
startTime: '13:00+01:00'
endDate: '2021-02-02'
endTime: '13:30+01:00'
youtubeEmbed: https://youtu.be/ZfxOFYXmiNw
locationURL:
  - https://doi.org/10.6084/m9.figshare.14160668.v1
  - https://youtu.be/ZfxOFYXmiNw
  - https://www.bilibili.com/video/BV1854y1h7d9
---

# nf-core/bytesize

Join us for the first webinar of a **weekly series** of short talks **“nf-core/bytesize”** starting Tuesday February 2, 2021.

The series will be very short lunchtime talks - approximately **15 minutes**, followed by questions / discussion.
They will typically be held every Tuesday at 13:00 CET.

## Bytesize 1: An Introduction to nf-core

This week, Phil Ewels ([@ewels](http://github.com/ewels/)) will present: _**An Introduction to nf-core.**_
He will introduce the new _nf-core/bytesize_ seminar series, give a short overview of the nf-core community
and describe some upcoming events.

The talk will be presented on Zoom and live-streamed on YouTube:

- Zoom: <https://zoom.us/j/95563936122>
- YouTube: <https://youtu.be/ZfxOFYXmiNw>

Questions/Ideas? Connect with nf-core on Slack: <https://nf-co.re/join/>

<details markdown="1"><summary>Video transcription</summary>

:::note
This text has been edited to make it more suitable for readers.
:::

Hello, thank you very much for the introduction.

It is great to get these talks started and speak to everybody today.
Apologies for the slightly last minute kind of announcement of these talks, they cracked up on us slightly, but it's great to see so many people joining today.

Hopefully this will be the start of a really exciting kind of series for us. The talks are deliberately really short and hopefully very focused. Today is a bit of a special one because it's the first one so I'm going to give a very short introduction to <b>nf-core</b> to the project for newcomers and then I am going to talk about the seminar series itself.

And also talk a little bit about the upcoming hackathon. If you have any questions then shout out or Renuka will curate those and we'll take them at the end so for any of you not familiar with nf- core

<b>nf-core</b> is a community effort to collect a curated set of analysis pipelines using Nextflow, the focus being very much on the community, really high quality workflows and working together to make that happen.

Nextflow is a workflow manager, a pipelining tool and the beauty of Nextflow is that it works on almost any computational infrastructure. The pipeline code is portable across all these different systems so you write one pipeline and it can run anywhere - on a cloud, or locally on your HPC.

It is this functionality which really allows a community like <b>nf-core</b> to work where we can build a single work together on a single pipeline and run it in very heterogeneous compute environments.

It is something which is extremely powerful as well, as Nextflow can also handle the different software container systems so you don't have to go and install hundreds of different bioinformatics tools every time you want to run a pipeline.

As long as you are set up with either Docker, Singularity or Conda or in fact also now Shifter and Podman and other containerised tools. All of the <b>nf-core</b> pipelines should work basically out of the box and each pipeline all comes with its own set of software dependencies already pre-packaged for you.

So that is Nextflow. <b>nf-core</b> is kind of a distinct project but of course we are closely tied to Nextflow. What we bring is a set of guidelines to harmonise all the pipelines within the projects. We are a large group of people now, and so we all need to kind of play by the same set of rules and the guidelines.

We provide a lot of helper tools which give additional functionality, both if you are running pipelines or workflow.

That helps to list what is available and download those pipelines so that you can run them on your offline cluster and kind of streamline some of these common workflows. But also for pipeline developers as well: to run automated tests or to build a pipeline from our template.

Finally what many people come to us for, is the resource that the pipelines themselves and so we have a large number of pipelines which are ready to go, some in development, still some some quite mature.

All of them hopefully useful and ready to go, the pipelines are great, the tools are fantastic but none of that happens without the community and that's really the heart of <b>nf-core</b>. What started off as kind of a small unofficial project mostly in Europe. This has really now spread to be fairly global and hopefully as inclusive as possible community.

As Renuka mentioned we are getting people from all over the world contacting us now and getting involved. This comes with its own challenges especially time-zones and things but it is great to see so many of you here.

The heart of this community and the communication is mostly Slack, so if you are not already, then hop onto slack and then you'll find channels for pipelines and also communicate and discuss all the different topics around <b>nf-core</b> and it's a really fun place to be. If you want more details on <b>nf-core</b> itself you can check out the publication on the [website](https://nf-co.re/publications).

What is this Bytesize? Bytesize is a new seminar series that we are starting off and nf-core has grown and more people were involved.

A lot of people have been wanting to get involved and get up to speed with how things work.
And there has been an increasing demand for training materials and documentation.
This hopefully is kind of another arrow in our quiver for that.

It is like a bit more material that we can start to generate training materials basically for newcomers. We'll be doing these as live seminars and that's partly to spread them out over time so that they're not too big an ask of any one person at any given time.

We will be spreading them out between speakers so it will not just be me talking every single week and also we have a question and answer session. Afterwards as well as that all of these talks are going out live to youtube where they will also be archived and curated into playlists so that they are easily accessible.

We are also thinking that we are going to put these in a prominent space on the <b>nf-core</b> website. And hopefully it might even start to associate the talks with videos there along with some written material where appropriate.

So this will be a kind of good way to take out specific focused topics, talk about them in detail and then if you're following along with the series you can learn as you go along. But also we can refer back to these talks if people ask us questions which are related and I have to thank Renuka for the excellent name suggestion and tagline which I think is very catchy.

We are aiming to try and do just 15 minutes of actual talking and it is because it is a kind of lunch-time seminar so we do not want to take too much time out of your day

We are going to try and do every Tuesday lunchtime in Europe, hopefully this time will work for other time zones as well.

We are going to try and keep each talk as focused as possible on specific topics so not too broad as many of our previous training talks have been.

We have been kind of the core team; the outreach team has been working quite hard the last few weeks and trying to think up different things. Our initial focus is going to be trying to cover some of the common training topics and topics which will be important for the hackathon.

This is really to try and bring everybody up to the same starting point so that when we do the hackathon you already know how to get started. This is an important note on the upcoming hackathon - previously we have had these meetings as a combination of training and actual writing code. But this time we're going to try and pull the training out a little bit, and spread that out over this seminar series instead and focus specifically on really like writing for the hackathon and make it a bit more a bit more intense.

So this seminar series is a good way to keep on track of everything and make sure that when you hit the hackathon you hit the ground running. We've got a lot of other ideas in addition to these listed here and this might change so we're going to be announcing these as we did today's talk. Hopefully with a bit more warning through the nf-core website.

We have the ‘[Events](https://nf-co.re/events)’ page, it was just suggested a moment ago that we make this a bit prominent. It is now listed as a main button in the top navigation which you can find easily that will list all the talks, when they are and what they are about.

We have a dedicated Slack channel on the nf-core Slack and you can pop in there to ask questions about any of the talks. We will try and scoop them up and of course we will be announcing them via our Twitter feed as well and we are really keen to hear what you'd like to know about.

So far the topics that we have come up with, the things that are commonly asked about within the nf-core Slack which are common topics of confusion or things which we think will help people. But of course we are going to want input from the community.

So jump onto that Slack channel and fire away with any ideas of things you would like to see short talks on and I will put them into the list and see if we can cover them. If you would like to give any talks we are are looking for more speakers always, so if you are if you are interested in giving a talk on any topic then please wave your hand, okay?

So that is the talk series. Also announced today is the nf-hackathon and this is our next big event as <b>nf-core</b> has grown over the past couple of years. We started off with really very small very informal hackathons and then they have kind of become more formalized, bigger, and more inclusive as we have gone along.

I am hoping that this next one will be a natural growth.

They have actually become one of my favorite parts of working with nf-core. Getting everyone together either in person or virtually and kind of having a real focused time to work together.
It is really evident then the community aspect because everyone is working together and chatting and you kind of see everyone's contributions.

Because of the ongoing pandemic the next hackathon will be hosted virtually, so it will all be online. It is going to be towards the end of March the 22nd - 24th.

We have done it at the start of the week with the intention that for people who are working on a project they will be kind of a tail off at the end of the week. So an unofficial hackathon can continue afterwards if you would like to.

Lots of good reasons to attend - of course we will all be there so you can come and hang out with us. This year Harshil suggested that we try and push this live pair programming which is becoming more and more possible with different code editors like Atom and Vscode.

So if we have two or more people working on the same chunk of code or the same idea, we can all hop into the same session and basically co-edit the same file live.

Even if it is just one person typing it is a nice way to actually kind of feel like we are sitting next to one another even though we might be on different continents. There is going to be quite a lot of focus on DSL2 (domain-specific language) too - the next iteration of the Nextflow syntax.

So if you have been hearing about this for a while but have never quite taken the plunge the hackathon is a great place to do. Because we are going to be doing lots of focused work on it, and leading up to that with some training materials. So it is a really good way to fully immerse yourself in this new world.

Speaking for myself as well there is also important stuff like Maxime has already started working on the collaborative playlist on spotify. So we can all get a groove on and of course um we are hoping to send out like we did for the last online hackathon.

If you sign up in time we are hoping to send out some little goodie bags with some <b>nf-core</b> stash, so that is reason alone. To be honest i'm not sure if it is 100% confirmed yet, but maybe a little something like this might be coming through your letterbox - some very fancy socks. Again thank you Renuka for this as we will see. No promises!

Just like the Bytesize we would love your help with the hackathon as well. The organization is a big job and it started already but there are some key things that we would really like extra people on board for. As this is online and we are becoming more global, it would be great if we could try to spread the hackathon out across different time zones and get more people outside of Europe involved - in the Americas or in Asia.

One of the ideas we had for that is just to have a handful of people who are in those different time zones or potentially willing to stay up late, so we know for certain that there is going to be at least one person around at any given time. Then we can try and keep the hackathon kind of rolling 24 hours which should be brilliant. Then if new people are hopping on and they do not know where to start or where to look, at least they know there is going to be someone responding to the messages and Slack.

If you would like to volunteer for that and if you're in a funky time zone that would be fantastic. We've also got the different hackathon projects, so we are looking for people to take the lead on those. For a similar reason, people coming in and wanting to know where to start or who to talk to, get a feel for which of the subtasks and what they could work on.

We are also looking for ideas for what people should work on different tasks and all of this is being organized through Slack.

So hop onto the new hackathon March 2021 Slack channel and either say that you would like to be involved or suggest any ideas. I should have put the link in here but there is also a webpage on the nf-core website for this.

Maybe I can share my screen window quickly just so you can see what I am talking about. Here you can see I am on the <b>nf-core</b> website. If I hop onto the '[Events](https://nf-co.re/events)' page , you can see that we have a page about today's talk. This one right now with links and everything.

There is also a page about the upcoming hackathon so if you click here, you will find the link to register. This is non-binding it does not cost anything to join a hackathon; it just helps us know roughly who is going to be there; also which projects people are interested in/working on, and if you want any of that cool stash then you have to put your address in so we can post it to you.

And then there are details about when it is going to be and how it is going to work. And this page will grow as the hackathon gets closer. You can also see some of the ideas that we are working on for different kinds of tasks to work on the project.

Great, that is my 15 minutes! So as I finish every talk, come and join our community if you are not already. we operate on Github, Slack, Twitter, Youtube and you can find all the details on how to join all of those on our website.

I did not introduce myself very much at the start but my name is Phil. I work at SciLifelab and I live in Sweden as well at the NGI. This lovely view I took with my new toy - it is about 70 meters above where I am sitting right now, when the sun was shining a couple of days ago

With that I am happy to take any questions. Please pop them into the Zoom chat or Slack, and Renuka can relay to me.

### Q&A

Renuka:
Thanks very much Phil and thanks everyone for listening. We have a question in the chat. Vidya wants to know whether there will be any talks on how to run <b>nf-core</b> on a grid system, like running on Slurm.

Phil:
Yeah I mean that is quite a specific question. Not so much nf-core, but that is more of a Nextflow issue we are going to talk about. One of the first talks will be about how to set up <b>nf-core</b> configs.
So we have a few different ways to bring in Nextflow configuration into the <b>nf-core</b> pipelines and to make that easy for you.

We have plenty of people using <b>nf-core</b> pipelines on Slurm, myself included and so that would be a good place to start. There is a talk coming up about that and also head over to Slack, we have a channel called configs which is designed for questions like this.

Renuka:
Okay thank you and there do not seem to be any more questions in our chat here.
Yeah so that seems to be it basically.

Phil:
I am sorry I meant to say I totally forgot I need to put on the finish line as well.
We need to put a big hands up and thank you to the Chan Zuckerberg Initiative (CZI).
These hackathons and these Bytesize talks are supported through a grant that was received from the CZI, the EOS programme (Essential Open Source for Science) project.

So Renuka who is organizing all these behind the scenes, she is partly funded by this and also
a lot of the backend stuff is helped through that. So it would be very difficult to do this without their help and we are very grateful. Remember to put that logo into the next talk!

There is also a question now about whether there is an agenda somewhere for the <b>nf-core</b> talks.
Yes, there will be so this will come up on the nf website under the ‘Events’ page.

We are just finalizing who is going to give which talks on which days.

As soon as we have settled on those, we are going to start adding them to the <b>nf-core</b> website, so check them out. I think when they are added to the <b>nf-core</b> website they automatically pop into Slack in a channel called #events so if you join that channel, you will see them appear there.

And we are also going to announce them on #Bytesize as they are coming up.

</details>
