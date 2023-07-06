---
title: 'Bytesize: gitpod.io'
subtitle: Chris Wyatt - Seqera Labs
type: talk
start_date: '2022-06-14'
start_time: '13:00 CEST'
end_date: '2022-06-14'
end_time: '13:30 CEST'
youtube_embed: https://www.youtube.com/watch?v=kBoC6QBU-M0
location_url:
  - https://www.youtube.com/watch?v=kBoC6QBU-M0
  - https://doi.org/10.6084/m9.figshare.20109701.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: gitpod.io

This week, Chris Wyatt ([@chriswyatt1](https://github.com/chriswyatt1)) will show us how to use GitPod (<https://gitpod.io>) for testing, development and training. GitPod gives you a complete cloud-based environment to work in for free, with VSCode editor integration (in your browser or locally).

Chris will show us how it works, how you can use the environments on offer and how to use the `.gitpod.yml` configuration that now comes as part of the nf-core pipeline template.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1)
Okay, hi everyone welcome to today's bytesize talk. Today we have Chris Wyatt who's going to be telling us about Gitpod and how it can be used for testing, development and training. I'll let Chris introduce himself but I hope you enjoy the talk.

[0:34](https://www.youtube.com/watch?v=kBoC6QBU-M0t=34)
Hi everyone, my name is Chris, I work as a postdoc at University College London and I also help with the training materials at Seqera Labs and this is particularly where Gitpod, basically a miniature computer in the cloud that we can use for development, is really handy. specifically for training but it is also used really frequently within the nf-core community to test code, make changes and use all the Git features of GitHub within this really nice environment.

[1:07](https://www.youtube.com/watch?v=kBoC6QBU-M0t=67)
What is Gitpod? It's an open source developer platform that's using VSCode and it can spin up these miniature development environments in seconds from a Git repository with all the code you need so that you can actually do all the testing or whatever you need to do in your particular piece of code. It's free to use and I'll get to more details of that later. All the pieces you need are actually just a Chromium based browser, an internet connection and a GitHub account. You will need to sign in using your GitHub account to get access to the Gitpod environment.

[1:43](https://www.youtube.com/watch?v=kBoC6QBU-M0t=103)
Just a bit over here what I'm going to talk about today, I'll explain why Gitpod is used and how you actually get started, which is a very simple setup, and then we'll get straight into running an nf-core test pipeline. Then we'll go into slightly deeper topics that are maybe more for developers of nf-core code, which is editing nf-core pipelines, dealing with GitHub and how to push them into branches for example, setting up these whole GitHub environments and live rendering of HTML and to do things with the website. There's also a really useful task that we use at nf-core.

[2:22](https://www.youtube.com/watch?v=kBoC6QBU-M0t=142)
Why are we using Gitpod? Well, it's really cool, it's the ethos of Nextflow, the fact that at nf-core that you have these repeatable and reproducible pieces of code but now we also have a way of testing them in the exact same environment between all the collaborators of a particular nf-core pipeline. It also helps to simplify a lot of the tasks that we're doing when we're coding with these particular pipelines and it's really fast. There are other different types of environments like this other than Gitpod, but Gitpod's super fast and normally you can just go straight from a GitHub repository to this kind of VSCode window where you have a terminal, all your code, and a really, really nice and neat way of dealing with all your work.

[3:06](https://www.youtube.com/watch?v=kBoC6QBU-M0t=186)
To get started we may as well just dive in, so that we can show exactly how this works. If you want to do this at another time, it might be difficult to do it as well at the same time but if you want to then go ahead. All you need is the browser and your GitHub credentials and to actually open any GitHub URL as a Gitpod environment. Then you just need to add this prefix that's shown here before the GitHub repository URL or you can install a Gitpod browser extension. I will put all of these in the notes afterwards. You can download this extension and it will basically add a Gitpod button, a green button, to every Git repository and all it's doing is just making this prefix before the URL, so that you can quickly access it with one click. Especially if you're already logged in with your GitHub credentials, then it will be one click and you'll get straight into the working repository. Finally, I just wanted to mention that this is actually being done through a file in the Git repository called `.gitpod.yml` file and this is the file that tells Gitpod what particular pieces of code we need and how to set up that environment so that the particular code that we're interested in works.

[4:29](https://www.youtube.com/watch?v=kBoC6QBU-M0t=269)
We're going to just go straight into the live demo. I'll explain exactly what it looks like when we go into this environment, how you can add extensions into the environment. You can add different add-ons. This is a picture from Gitpod. You can add Docker or you can add lots of different functionality within your containerized environments to test your code. I will explain how to use Git within this and also opening files like PDFs etc and live rendering for website development.

[5:02](https://www.youtube.com/watch?v=kBoC6QBU-M0t=302)
What we're going to do now is, we're just going to go straight into how to use Gitpod. Hopefully, everyone can see this. The example pipeline we're going to use is the nf-core test pipeline. I've already got the browse extension so I actually have this link already. I could also copy the link address, open up a new tab and just paste it in. But for speed, I'll just click on it. So when you log in, all you need to do... and at this point, it would ask you for your credentials. For me, it's not because I'm already logged in. This is going straight from that repository into Gitpod. On the other side you can see the explorer set says, these are all the files within your repository; you have a window here to look at visualization of code or PDFs etc and then you have a terminal at the bottom. You also should have extra buttons here which I can't see. Anyway, there's all the things you would normally want to do on your own computer, you can also do here.

[6:21](https://www.youtube.com/watch?v=kBoC6QBU-M0t=381)
Just to show you how this is actually working, the Gitpod.yaml file is not... Did I click the right one? I'll go back. This particular repository is in the master. For some pipelines (and this is a very good example), that you want to run on Gitpod, you need to make sure it's got this Gitpod.yaml file, because this is all the programs you need. If we click on that, it's showing that we have an image, so a docker file from dockerHub. The image is being pulled, that has all the code that we need and I will go a bit more into these Gitpod yamls a bit later.

[7:05](https://www.youtube.com/watch?v=kBoC6QBU-M0t=425)
We need to be in this development branch. I also mention at this point, if you're actually making changes to nf-core code then you'll be wanting to do this in the dev branch and not on the master branch. If we click it, we should get back to a Gitpod environment, that I was expecting. That should have all the code that we need to run this particular pipeline. Sometimes this can take a few seconds. Probably because it's a live demo, it'll be more than that, but if not we can move on to different things I can explain to you. It should pop up reasonably fast. It won't? No. Okay, there we go. Now we are actually in the correct repository that I wanted to be in, but just move that away, just to make it clearer. You can see there's a gitpod.yml file which I was expecting to see. This is the reason why we have the code brought into this environment. In the terminal, if we type in `nextflow` it finds the tool because it's already installed within this gitpod.yml file. It told us that this container contains all the things we need to run the particular piece of code.

[8:18](https://www.youtube.com/watch?v=kBoC6QBU-M0t=498)
The other thing that's quite nice here if you want to actually test this repository, there's obviously the readme. If we open this... This is obviously not rendered... There's a preview button, so we can easily click things to preview the render to look exactly the same as it would do on the web page. It also tracks, so if you're actually needing to change something in the readme you can add whatever you want and it will automatically change it in the preview on the right-hand side. But I'm not going to change this at the moment. If I move this across, if we look at the readme, it tells us how we can run this piece of code and I just want to quickly show you that it's quite easy to quickly run this particular piece of code.

[9:05](https://www.youtube.com/watch?v=kBoC6QBU-M0t=545)
We're going to run nextflow. Remember that this is actually in a dev branch. If we run it exactly the same as this line, we would actually pull the master branch, which is not actually what we need to run. So obviously with nextflow we can just run from the main.nf file. Then we can set up a profile. Here it explains how to set the profile as normal with nextflow code. In this particular case, we're going to use docker. The reason for this is that singularity is one of the few things I know that doesn't work in Gitpods. If you need to test something with singularity you'd have to do it outside. You wouldn't be able to do it within this environment. Then we could set an outdir. Everything you need to do here is very simple to do. If you want to add a new folder, we could call it results. This is the way you can add a folder manually, by clicking one of these buttons for new files on your folders and then we can just run the code... like this.

[10:11](https://www.youtube.com/watch?v=kBoC6QBU-M0t=611)
Now we're running Nextflow and this pipeline within Gitpod. I hope you can see it's quite nice that we could quickly get to run this pipeline, see how it works and begin to see the output. The resources on these machines aren't huge, there's only a few threads, so you can't do anything that's going to be using, like, eight threads for example. They only run for a maximum of half an hour if you have the free account, but we'll explain later, that there's an access where you can get a limited time if you need that for your particular project.

[10:46](https://www.youtube.com/watch?v=kBoC6QBU-M0t=646)
If we just see that the program is running and we'll see all the tasks connected are running through, we can then... I'll just wait for it to finish, shouldn't take too long... This pipeline, I should have mentioned, is actually taking FASTQ reads. It's running it through FASTQC to get quality control of those particular samples and it's using MultiQC to push them all together to make a really nice plot of all the different features and the quality scores of the FASTQ reads.

[11:20](https://www.youtube.com/watch?v=kBoC6QBU-M0t=680)
This is a quick way of running a pipeline and you can do this for any of the nf-core pipelines that have this gitpod.yaml. It's quite fun to do. It means you can just really test out a new pipeline and see exactly how it works and if it works with your data. You can drop and drag your files into this explorer area to bring files in. There's also other ways that are normally here. There'll be a file in store, I don't know why it's not here. Might be because of the screen. Anyway, you can drop and drag files in here. I'll just take that off just see if we can find it. Yes, here we go. Sorry if I have the actual full screen. You have options here that you can save files to your local machine and you can open files and do all the kinds of things you would expect to do on a normal computer.

[12:06](https://www.youtube.com/watch?v=kBoC6QBU-M0t=726)
The other thing, now it's finished. In our results file, we have our results and if we look in MultiQC, for example, we have the MultiQC report. Again, you can use these preview buttons, so that you can render these particular results into a nice file where you can see everything that's going on. Here is the results of that FASTQC result in MultiQC. I think this is quite cool to know that you can just quickly run these pipelines.

[12:36](https://www.youtube.com/watch?v=kBoC6QBU-M0t=756)
Within this particular environment, you might want to change something. So maybe you're particularly interested in changing a particular module. You can still change everything in a Gitpod environment and once you change something... I'll just add... I won't actually push this to the branch, because people won't be happy. But if I add on `cat` or whatever I've changed here. If I'm changing the code here you'll notice on the left-hand side, on the source control, it's noticed that one of the files within your repository has changed. If we click this button it said that we've changed the main.nf and you can now use the git functions within Gitpod to add that particular change and you can also then make a message for that particular change. Also, down the bottom, you can see that we have this development branch - on the bottom left-hand side. If we click this we can also choose to create a new branch. If you've made changes to the code this is a very quick way to create a new branch and show that we can push it back to the repository and then do a pull request.

[13:47](https://www.youtube.com/watch?v=kBoC6QBU-M0t=827)
I think that's everything I really wanted to show you in this very basic example. I'm going to quickly go to another example that I think is quite helpful to explore, which is actually on the website. If I go to the website repository we can also run it from here. I'll just show a few more features that I think are really interesting to see. We can actually live render the whole website and show that the whole website is changing and we can change the underlying markdown. This is quite helpful when we've been developing the website so hopefully, this will quickly load and I can show you this example.

[14:31](https://www.youtube.com/watch?v=kBoC6QBU-M0t=871)
I will show you this repository. We see the gitpod.yml file, this is what I showed you before. There are many things you can put in these yml files and I think the best way to learn about this is to go to the gitpod.io website or by looking at other people's gitpod.yml files. There you can see exactly what people are doing.

[14:56](https://www.youtube.com/watch?v=kBoC6QBU-M0t=896)
Before I showed you `image`. Image was pulling a docker container, but you also can initiate tasks. You can do lots of different commands related to docker. You can open ports to expose ports so that we can look at web browser and share a particular web browser HTML rendering. The other thing is the VSCode extensions. These extensions are all over here and you can make sure that each repository that you make into a Gitpod environment has all the different code linting services, adding different particular pieces of code that are very useful to have each time. You can add these extensions manually in this gitpod.yml.

[15:43](https://www.youtube.com/watch?v=kBoC6QBU-M0t=943)
You can also set all the GitHub settings that you would need to change. Just very quickly, because I know I'm probably running overtime, here in the remote browser we can set this to being a live preview. Okay, no you can't... Probably I have done something wrong. There we go! Within this repository, we're now making a live view of the actual repository that's changed and then if I wanted to go through the code into the markdown I can actually change the markdown. Like I showed you before you can actually change what's happening in this live preview of the nf-core website. I think this is super useful as well.

[16:30](https://www.youtube.com/watch?v=kBoC6QBU-M0t=990)
The last thing I wanted to show you was with the GitHub button here, you can look at all the open pull requests you can click on different ones that you're interested in trying to fix. We'll get rid of this for now. It should show you the pull request and give you all the information. If it loads... hopefully it will... Anyway, normally it should come up with all the things you can do on git basically you can do within this Gitpod environment, I think this is really handy. So there we go!

[17:04](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1024)
You can do everything that you want to do in git, but it's within this GitHub repository and I think this is really helpful. Finally, I think I showed everything I needed to show there. I just wanted to mention it here... it's only a few slides left... I wanted to mention that also these Gitpod environments are used for extensive training. Sateesh presented really well all the training materials that we have two weeks ago. All of these repositories are now running a Gitpod environment, with all the code that you need to run them. If you're new to Nextflow, this is a really cool way of getting into Gitpod and Nextflow, that you can have this amazing environment in the cloud with all the tools that you need.

[17:48](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1068)
The cost. The free one is really cool, it's only 50 hours per month. You have four parallel workspaces, but you have this 30-minute activity time. As a member of the nf-core GitHub account you can find a free professional open-source account. You can get this through the following link. If you click this, there's a form you need to fill in to be able to get a professional account. You can get this through the link. If you click this, there's a professional account for free, which is really nice of Gitpod. For open source developers. Just to say there's not very many threads and I think the maximum size of each repository is about 30 to 50 gigabytes, in case you're testing and you need more than that. But normally the resources are easily enough to do all these different testing things that we want to do.

[18:32](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1112)
Finally, if you need help, there's lots of people in nf-core that work with gitpod and have lots of knowledge of how to run all these things. How to make the yml's for example. There's also a section on the web page relating to Gitpod and there's really great resources and videos on gitpod.io itself. I think that's everything. Thanks for listening.

[18:56](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1136)
Thanks Chris, that was a really fantastic talk, and really exciting to hear a bit more about Gitpod again. Are there any questions? The options are to either put your questions into the chat and I can read them out, or I will just unmute the microphones so you should be able to speak up yourself.

[19:27](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1177)
(question) Security issues (from chat).

(speaker) Personally, I really don't know, I would have to check that up for you, I don't know about, what level of security they have.

(comment) I would say for all of the stuff to do with nf-core there's basically no security issues. It's public code running on a public server. If you want to do something secret on a private repo, then maybe consider that. But nf-core stuff, I don't think you need to worry. It's completely off your system, so it's more secure than if you're running it locally.

(speaker) True.

[20:11](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1211)
(question) what are the pros and cons of Gitpod compared to code spaces

(answer) This is a good point. We did play around with a few of these different things.Gitpod is super fast. I think it's powered by kubernetes and it's done in a much more efficient way. Iit means that everything's really quick and you have these initializations step. The reason it's quick is because it actually saves these layers of these Gitpod environments so that you don't have to keep installing everything at the same time. That's why it's like a second, or a few seconds normally, to get to these things. It's mainly for speed and also they're a really great company that seem to have a really active community of people that help you with your code. That's the main reason.

[20:52](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1252)
(question) The next question is, is there a way to put a link on a document to open it in Gitpod, that automatically takes the current branch into account?

(answer) I assume so. I'd have to check that again. I'm not sure. I'd have to try it but yes, I think it must be possible. I just didn't didn't share that. I was meant to share the links, I've just done it quickly.

(comment) You can just prefix the file link with the gitpod.io. I think it works.

(speaker) thanks Mahesh

(comment) It takes into account your entire current URL. If the current URL you're looking at is on a branch, you've got to open Gitpod in that branch. And the same with the little green Gitpod button. Also most of the repos we have - a lot of them now - when you open a pull request Gitpod has an integration, which will automatically add a little button into the bottom of the first description of that pull request. When you're browsing pull request you'll also see a button saying "open in Gitpod" and that will open you directly into that pull request.

[22:01](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1288)
(question) There's another question here. Can you kindly share the links e.g the slack channel link here?

(answer) Yes we will make sure the links are at least on slack on one or two channels to find later.

[22:19](https://www.youtube.com/watch?v=kBoC6QBU-M0t=1339)
(host) Unless there are any more questions that pop up quickly, I want to thank Chris again for his fabulous talk. You can always catch up on this talk again on youtube. We'll post this very shortly after the talk. We'd also like to thank the Chan Zuckerberg foundation for supporting these talks. Thanks very much.

(speaker) Thanks guys!

</details>
