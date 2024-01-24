---
title: 'Bytesize: transcripts of bytesize talks'
subtitle: Franziska Bonath - NGI, Stockholm
type: talk
startDate: '2023-01-31'
startTime: '13:00+01:00'
endDate: '2023-01-31'
endTime: '13:30+01:00'
youtubeEmbed: https://www.youtube.com/watch?v=amwwmFMwOYw
locationURL:
  - https://doi.org/10.6084/m9.figshare.21995243.v1
  - https://www.youtube.com/watch?v=amwwmFMwOYw
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: transcripts of bytesize talks

This week, Franziska Bonath ([@FranBonath](https://github.com/FranBonath)) will talk about her work to generate transcripts of bytesize talks and what these might be used for in the future.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=amwwmFMwOYw&t=1)
(host) Hi, Maxime here. First of all, I'd like to thank the Chan Zuckerberg Initiative to help us doing these bytesize talks. And today, Franziska Bonath will present us how the transcript of the bytesize talks happen. It's a very meta bytesize talk today. And as usual, please use Slack for your questions. Now, it's up to you, Fran.

[0:28](https://www.youtube.com/watch?v=amwwmFMwOYw&t=28)
Okay, thank you. Welcome, everyone. I'm talking about bytesize talk transcripts. Just very briefly, what we're going to do today. I will handle the question of all questions, why to transcribe bytesize talks at all, and then briefly go in how we did it and why, what we're going to do in the future.

[0:52](https://www.youtube.com/watch?v=amwwmFMwOYw&t=52)
Why, why are we going through all this pain? And the big answer is that we want to be more inclusive. This is one of the reasons why we got funding from the Chan Zuckerberg Initiative. But of course, we also have a desire to do this. Not everyone is able to hear things. If you rely on the transcripts that are automatically done by YouTube, for example, it can be very difficult to get the gist of what the talk is about. Also, even if you hear perfectly, not everyone will be able to understand English well enough to figure out what the talk is about. In addition, we have the speakers from all over the world. There might be accents that are a bit more difficult to follow. And so, having a really good transcript will help understand these talks a lot better.

[1:54](https://www.youtube.com/watch?v=amwwmFMwOYw&t=114)
There's other reasons. One is, of course, to improve the subtitles for YouTube. But also, if you have the transcript in itself without the video, you should be able to understand it. And it will be a resource for understanding of details that are maybe not in the slides. There will be, hopefully at least, the correct names of all the tools that are used. You can look that up, and then it will be easier to search for that online. But also, once you have a text, there's a lot of things you can do with that text. You can translate the text, you can put it into some AI-based thing and have it give you a summary of the text. There's a lot of things that we might start to think of in the future, and it will be text-based. And the better the information is that you give in, the better it is what you're going to get out.

[2:55](https://www.youtube.com/watch?v=amwwmFMwOYw&t=175)
Where can I find these transcripts? It's, at the moment, a bit difficult. I admit that. I'm going to quickly show you. What you have to do at the moment is, you have to go to... [...] If you're on the website, you go to events, and then you can search for only bytesize here. And this will be the upcoming ones. But if you go to the past ones, for example, let's go to taxprofiler, and you scroll down. What you will find here is the embedded YouTube video. And at the bottom, you will have the transcripts. You can go directly to one of those. It will show up there. And this is, at the moment, the only way how you will get the transcripts for any talk. But it will be uploaded to YouTube eventually. Then I go back to my slides.

[4:13](https://www.youtube.com/watch?v=amwwmFMwOYw&t=253)
How did we do this? We did try to use the automated transcripts on YouTube first. It is horrible. Basically, what happens is that you have a lot of these oohs, and aahs, and ohms that are not removed at all. Also, you will have no punctuation whatsoever. You have to add the capitalization after every full stop that you have in your transcript. It takes forever. It probably would have been quicker to just write it while you hear it. That did not work. And that means in comes a new tool, which I'm forever grateful to Matthias Zepper, who introduced me to it. It's called Whisper. And at the moment that I started this transcript, Whisper was only available as a tool as is. But from now on, you can also have a Nextflow pipeline for Whisper. You can find it under this link. And Whisper helped with a lot. It does add punctuation. It does surprisingly recognize a lot of the tools that we're using. And it removes all the emms. It removes a lot of the double mentions. If you're talking normally, often you stop for thinking about something and then you repeat what you have just said before. And so, these double mentions, they get edited out automatically, which is super nice. I can only recommend Whisper if you ever do transcripts of any video yourself.

[5:51](https://www.youtube.com/watch?v=amwwmFMwOYw&t=351)
But even though Whisper is great, it is not perfect. I don't think any automated talk transcript ever will be perfect. The main things that we have to do is add timestamps so we have nice sections that belong together. But of course, also names, specifically names of people, but also of tools often get not identified correctly. You have to check and edit those. Specialized terminology is also not recognized because they are probably not in the library of Whisper. And also sometimes sentences are super long. It might be ellipsis or that someone had a thought, stopped in the thought and continued afterwards, which is totally fine if you're just listening to a person. But if you want to just read it, it's very difficult to understand. These kinds of things we have to manually change afterwards.

[6:55](https://www.youtube.com/watch?v=amwwmFMwOYw&t=415)
To give you a kind of an idea. Our most favorite words that are part of pretty much every bytesize talk, nf-core and Nextflow are very commonly misspelled. Nf-core very typically gets misspelled to NFL and NF4. I don't exactly know why, but in every third or so transcript I read those. And of course, you also have just some misspelling of nf-core itself. Sometimes it does pick it up and in very, very rare cases, it will also type it correctly. Nextflow, it also has diverse ways of how it can be written. In the latest one, I had it transcribed to "next floor". But then of course, there's just some random things that don't repeat. Like, "elusion" will be transcribed to "illusion", "iterations" to "situations". One of my favorites was "bioinformaticians" to "by partitions". Surprisingly, bioinformaticians, which is not that uncommon a word, I would say, gets transcribed a lot wrong. And you can imagine that if you have ribosomal RNA mistranscribed to rivals of RNA, the sentence will not make any sense. The handy overall summary can also become a handy oral summary, which would make sense, but which would change the meaning a bit. And just one other example, if you have a sentence like, "these processes take a sort of BAM from the samples", if you just read the sentence, I would have not guessed specifically what this would mean. Once I listened to the transcript, it turned out that it means "these processes take in a sorted BAM from SAMtools". This, I think, shows very clearly that manual work is necessary and that it's worth going through this and make these changes and not just rely on an automated transcript.

[9:07](https://www.youtube.com/watch?v=amwwmFMwOYw&t=547)
Now we're done, right? We are up to date. Everything's fine. Not quite, obviously. We have to add these transcripts to the subtitles on YouTube, which will happen in the not too far future, I hope. And also, what we want to try and see if we can do translations of these YouTube transcripts that we generate now, to have them in different languages, which would be super nice. Of course, bytesize talks are not finished yet. In fact, this very bytesize talk is going to be transcribed. We have this kind of Inception way where a bytesize talk that talks about bytesize talk transcripts is going to be transcribed. Anyway, this was all. I would like to thank Matthias for his enormously helpful tip for Whisper. And also, he was writing a container, I think, for Whisper. Marcel and Christopher, who had to approve all my pull requests for the transcripts. Of course, all the other reviewers, specifically the speakers that went through this horrendous task of reading their own talks. I'm not looking forward to this. Thank you very much. Now I'm open to any questions. Of course, there's no repository nf-core pipeline. I just took the... Anyway, thank you very much, everyone. Off to Maxime.

[10:53](https://www.youtube.com/watch?v=amwwmFMwOYw&t=653)
(host) Good. That was brilliant. Thank you very much. I will try to allow everyone to unmute themselves if you have questions. We haven't done that in a while. Where is this?

(speaker) Yeah, this is what you get when you use the template. Unnecessary things get included in the talk.

(host) Is there any question, actually, like, oh, yes.

(question) Jasmin is asking, how was the transcript added to YouTube? Will they be visible as normal subtitles?

(answer) So, yeah, I did look a bit into that. You can add your own subtitles in YouTube if you are the owner of the YouTube channel. As nf-core, I can add subtitles, and it will be one of the different subtitles that you can choose from. I think it's going to be called... No, I don't recall how it's called. But I think it will be the default option as subtitles.

[11:56](https://www.youtube.com/watch?v=amwwmFMwOYw&t=716)
(host) Okay. Do we have one last question, or are we good for today? I think we are good for today. Thank you again, Fran, for this presentation. Definitely that was a question I had, how everything was happening and stuff. Thank you very much for inviting me into all that. And see you soon. Thank you.

</details>
