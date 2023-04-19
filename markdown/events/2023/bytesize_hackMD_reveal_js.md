---
title: 'Bytesize: HackMD and reveal.js'
subtitle: Maxime Garcia, Seqera Labs
type: talk
start_date: '2023-04-18'
start_time: '13:00 CEST'
end_date: '2023-04-18'
end_time: '13:30 CEST'
location_url:
  - https://kth-se.zoom.us/j/68390542812
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/HackMD and reveal.js

HackMD.io is an open source, collaborative Markdown editor, which allows people to share, comment, and collaborate on documents. At nf-core we use HackMD extensively during hackathons and for notekeeping.
This week, Maxime Garcia ([@maxulysse](https://github.com/maxulysse)) will show us how to create Presentation Slides on HackMD using the `reveal.js` integration.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://kth-se.zoom.us/j/68390542812&t=1)
(host) Hello, everyone, and welcome to being back at the bytesize talks. Since we had a bit of a break over Easter. Today we start the series with Maxime from Seqera Labs, and he is going to talk about how to use HackMD to present slides. Off to you.

[0:24](https://kth-se.zoom.us/j/68390542812&t=24)
(speaker) Thank you, Fran, for the introduction. Hello, everyone. I'm Maxime Garcia working at Seqera Labs, and I'm going to present one of the tools that we use quite often during the hackathon and everything. Which is HackMD and reveal.js. Let's start with the presentation mode, because I think this is the most common stuff. The usual disclaimer, I'm mainly covering my own usage within nf-core or what I usually do on the side, but you can do much more. Don't hesitate to investigate and explore, I think it's fun. Yes, this is messy, but it's fine.

[1:09](https://kth-se.zoom.us/j/68390542812&t=69)
First, Markdown. Because all of the stuff that we're doing is in Markdown. What is Markdown? It's a lightweight markup language. It means that it's using small tags to do stuff, but most of the time, it uses just tiny symbols as the tags. The key point is readability. Because if you can read it... this is markdown... it's understandable, at the same time, what you're seeing here is what you see there, and you can understand what you see and what you read. For me, that's one of the key points in Markdown, it's readability. If you compare that to LaTeX or some other language, which are less readable, Markdown is super high in readability. It's super easy to convert into HTML, and then PDF, and it's widely used in all of the nf-core documentation and the nf-core website. You've noticed some tiny issues around the site, but I will come to that later. Don't worry about it.

[2:18](https://kth-se.zoom.us/j/68390542812&t=138)
Quick links for Markdown, I will of course share the slides later. It will be super easy to follow the slide to follow the links and everything. Quick links, I think I have one first quick reference to the Markdown syntax, which is fairly simple, a more complete Markdown cheat sheet, which goes way more into detail. Then some more documentation for the GitHub-flavored Markdown, which is an extension of the regular Markdown.

[2:58](https://kth-se.zoom.us/j/68390542812&t=178)
What is HackMD? HackMD is a tool for real-time collaboration of Markdown documents. It's widely used during hackathons and nf-core bytesize presentations. It means I can edit a document. Just like that. Real-time edit. If someone else has the same URL as I, we can all edit the same document, exactly what we do with Google docs. It's all in Markdown, so it's super easy to do that. HackMD has the possibility to use reveal.js, which is an HTML presentation framework, which we also use widely during the hackathon and that I use a lot for the bytesize presentation. Reveal.js is another tool and HackMD made it possible to use reveal.js, in it's presentation mode.

[4:00](https://kth-se.zoom.us/j/68390542812&t=240)
What you do when you have your presentation, what you can do, you can share in slide mode. If you share in slide mode, then you can have your presentation directly as a slides, and that's all. To pass from one slide to another, once you're in slide mode, you either use the arrow on your keyboard - I'm pretty sure it works also with arrows here, or yes, it does - or you can also use the space bar. Reveal.js allows you to do sections and subsections. You can do a lot of stuff within reveal.js. You can use fragments, so if you want to have a multiple steps in the slide that get revealed one after the other, or stuff like that. But usually what I do is super simple slides.

[4:53](https://kth-se.zoom.us/j/68390542812&t=293)
How to actually use reveal.js to present within Markdown. First it's super simple. I linked how to create a slide deck, which is the HackMD documentation for it. The most important part is to follow syntax, and the most important part in following syntax is to separate the slides. You need to have one empty line, three dashes, and another empty line. That's all. Then I usually use one or two hashes, so h1 or h2 or h3 for details, subtitles. That's all, you do as you want. Sometimes you might need small text, I usually use HTML tags. Small text, or big if you need big, but usually I just need to make things smaller. What I do also, sometimes if I need some really small text, I use also font. I think there is a font tag that you can specify the size of the font. That's super simple.

[6:06](https://kth-se.zoom.us/j/68390542812&t=366)
If you want to include a picture, if it's available online, you can link it as you usually would link any picture in Markdown. You just need to know the syntax for it. You can use HackMD directly to upload your picture. If I click there, I go to "insert image", go into my download, I can insert one document, and here it is, my document should be inserted, and here it is, I have it in my image line. As you can see, it's super easy to upload a file from your computer directly to HackMD, and this is the syntax to have a picture in Markdown.

[6:55](https://kth-se.zoom.us/j/68390542812&t=415)
In HackMD, we can also use fontawesome, which is what we use a lot in nf-core to have some... I know what is causing that, so I will explain that later, sorry about that... fontawesome is resources that allows you to use simple icons for different stuff. You can have, here it's for server, but you can have a thing that was also, GitHub should be something well known, Facebook should work as well, and so you have a lot of different things: cloud, cloud should work. If you know it, it's super easy to use. What we can do, but that is regular, so yes, simple. Let me put back what I did before, okay.

[7:57](https://kth-se.zoom.us/j/68390542812&t=477)
Background, this is super useful whenever you're doing a presentation, you might want to change the background of your presentation. What I usually do there is that I use this dot slide command, and I specify which is my data background. We can add some opacity if we want by using a data background opacity. I use here 0.5, but you can have more or less, so 0.8, 0.1, so you can really play with that. Let's put it back to 5.

[8:34](https://kth-se.zoom.us/j/68390542812&t=514)
An important thing also we can do with reveal.js is export as PDF, so this we can do only when we're in slide mode. Let me show that to you.. When we're in slide mode, if you scroll down, we can see some links. We can see that I made some changes a few seconds ago, I can see that nf-core is owning this note, because of course I made this note in the nf-core organization. I can edit back the note or I can also print the note. If I print it, it's going to try to print everything as it is, and what we want actually is to print everything as a PDF, so I'm printing to file as a PDF, and okay, and it's printing. I should have my file saved on my side, and it's efficient because you really get one page per slide. I've noticed from time to time some issues. Mainly if you have pictures, like here, that are way too big. In that case I would recommend to go back and forth in between from what you have on your screen and what you have in your PDF, and try to scale down your picture or have your picture on the next page. But apart from that it should work fairly well. Most of the issues that can happen in this case is because you're not following the proper syntax.

[10:16](https://kth-se.zoom.us/j/68390542812&t=616)
What we can do more with HackMD, you can include your own CSS slide, so your own CSS style, this is what I did here, and I'm guessing this is what is causing me the issues that I have there. I just have been trying that out yesterday, and yesterday I didn't have this kind of issue, so I don't know where it's coming from. I will have to look more into that. What you're doing is just saying that I want HackMD to include this file from HackMD, and this part is exactly correspond to a document there, so this is fairly good. You can also use that to include a simple slide deck. For example what I'm doing here in this presentation, I'm including the last slide deck, and what I'm doing if I go on the page for this slide, I can see directly what is linked here. I can see here this is my own reveal.js style, and I have all my style into HTML tags, and all of the CSS correspond to this slide, which could help you. What I really like about doing that is that you can add your own custom stuff, and what I like about doing that is adding a multi-column possibility, and that helped me using "div". I'm using divs within a div multi-column to present a different... here I'm just presenting three lists into three columns, and that's all. You could really do whatever you want.

[12:08](https://kth-se.zoom.us/j/68390542812&t=728)
Then I have some more tips. When in doubt, if you have any issue, don't hesitate to add more break lines. I've noticed that it's something that works a lot in my case, because I noticed that reveal.js and HackMD, the combination between the two is very, very dependent on the break line and the syntax. It expects you to follow the proper syntax, and sometimes it's not working very well. If something is not working reload or relaunch. If you can, you can also lint markdown if you're thinking that something is not working. Printing works much better when you're using Chrome, and also following the proper syntax. The less HTML, the better as well. I think the div works well, but apart from that, I'm not sure. I think I noticed some issue with tables and stuff, it's something that you really need to look for. What you can do in reveal.js, you can set up your username. Here in my case, my username is @maxulysse, and that links me, directs me, to myself. Otherwise, if you haven't set up your own username, you will have a something which is more or less un-understandable. As I explained, you can also link directly another slide, which is what I did just here. I think I'm good for questions, so I can open everything.

[13:53](https://kth-se.zoom.us/j/68390542812&t=833)
(host) Thank you. Now anyone who wants to ask a question can do so. I enabled you to unmute yourself. Are there any questions from the audience?

(question) I had a question. Is it possible to add video files to reveal.js syntax?

(answer) Yes, definitely. reveal.js is just HTML that is presented within JS, so you can link YouTube stuff. I don't do that often because most of the time I'm assuming a presentation should be shared easily in other settings, so I usually don't include videos, but you can do that. I think I did that in a previous presentation. I think other people have done that. Definitely, I think if you have a look at all of our presentation on nf-core, you will find some that link videos.

(question cont.) Okay, thanks.

[14:59](https://kth-se.zoom.us/j/68390542812&t=899)
(host) Are there any more questions?

(speaker) All right, this is annoying. I will not include the CSS anymore in the stuff. That was super convenient because you don't need to have your own style, but yes it's... or maybe I will try to figure out what is unreliable there.

[15:19](https://kth-se.zoom.us/j/68390542812&t=919)
(question) Maybe a question from my side or clarification. This is a free software, right?

(answer) Yes, definitely this is a free tool. I think we have something special, because we have an organization, so we might have something there. Oh yes, we have a team plan. I have no idea what it means for nf-core, but we have a team plan there.

(question continued) But if you would want to use it privately or...?

(answer continued) Oh yes, you can use it privately. Like I use it privately as well where I'm making my own list of books to read or that I will share with friends. Yes, you can make your own stuff exactly, I have a couple of presentation as well and yes you can have your own workspace and everything. You can have a private notes and the public notes. Yes, you can do some stuff.

[16:19](https://kth-se.zoom.us/j/68390542812&t=979)
(host) Perfect, thank you. Ah and there is also a link from James which explains the difference, I guess, between free and paid version of HackMD.

(speaker) That's good because I guess I had definitely no idea what is the difference between the two of them.

(host) Okay cool, thank you very much. Are there any more questions? It doesn't seem so. Then I would like to thank Maxime again for the nice talk today and all of you for listening and of course as usual the Chan Zuckerberg Initiative for funding our bytesize talks. Thank you very much.

(speaker) Thank you.

</details>
