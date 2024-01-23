---
title: 'Bytesize: nf-core/pre-commit'
subtitle: Matthias Hörtenhuber, Scilifelab Data Centre, Sweden
type: talk
startDate: '2023-02-07'
startTime: '13:00+01:00'
endDate: '2023-02-07'
endTime: '13:30+01:00'
youtube_embed: https://www.youtube.com/watch?v=08d6zv6zvdM
locationURL:
  - https://doi.org/10.6084/m9.figshare.22047485.v1
  - https://www.youtube.com/watch?v=08d6zv6zvdM
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/pre-commit

This week, Matthias Hörtenhuber ([@mashehu](https://github.com/mashehu)) is going to explain the use the newly added pre-commit tool added to nf-core/tools. Pre-commits were developed to inspect the snapshot that is about to be committed and helps to check formatting etc. before adding the code to the repository.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=08d6zv6zvdM&t=1)
(host) Hello everyone to today's Bitesize talk. I'm Franziska Bonath. I'm the host today. With me is Matthias Hörtenhuber. He is going to talk about pre-commit.

[0:11](https://www.youtube.com/watch?v=08d6zv6zvdM&t=11)
Yes. Hi. It will be a short one. Let's get going with it. The title is "pre-commit, hooked on a commit". It's about how we're trying to solve this problem a lot of us run into: You commit your code, then GitHub CI runs prettier, it finds a mistake and asks very passive aggressively, if you maybe forgot to run prettier. First solution is, of course, to ask nf-core-bot to please fix the linting. I did a quick search. Actually, 260 PRs have at least one comment, where somebody asked the bot to fix the linting. That's nice that people are using the bot. But actually, maybe we should try to save the ice bears and give nf-core-bot sometimes a break, because these things we can do differently.
What if you don't need to run the bot, because you already run prettier when you commit your code, or actually before you commit your code? We use this tool called pre-commit, which runs prettier when you hit `git commit`. It runs prettier, checks if the code is fine and makes changes on it.

[1:43](https://www.youtube.com/watch?v=08d6zv6zvdM&t=103)
How do you set pre-commit up? Well, the good news is it comes already pre-installed with the nf-core tools as a dependency. Then you need to have a `pre-commit-config.yaml` file. For example, this one here for prettier. Also, more good news, in the next tools release, this will be part of every pipeline template. Also in the modules repository, we have that set up. It doesn't change anything for you that we have it there. But if you then also run in your repository `pre-commit install`, it actually sets up this git hook. Whenever you hit commit, prettier is run beforehand and doesn't allow you to commit until you fixed these changes.

[2:39](https://www.youtube.com/watch?v=08d6zv6zvdM&t=159)
How does it look like? I made here a short example where I just added a line in the README file. As you can see, it's just below a heading, so Prettier will not like it. I run commit. Actually, prettier was run and fixed it. But that's the important thing, it didn't commit it yet. It's changed but not added. I actually need to run commit a second time. That's something I sometimes forget. You always need to run git commit twice if there is something wrong. If nothing is wrong, your prettier passes, then the commit runs through. This is prettier, which we use for markdown files and similar files.

[3:35](https://www.youtube.com/watch?v=08d6zv6zvdM&t=215)
One of the nice things with pre-commit is that you can use every code linting tool you want. Also, it doesn't matter if it's in a different language. Like you see with prettier, which is actually an NPM tool. We don't need to have node installed to run this version of prettier with pre-commit. It just comes through the mirror there. But other tools like Python-based tools come directly from the tools itself. If we set up this config, like we have for the tools repository, for example, it automatically checks the Python files with Black and isort. In this example, I added again to the README file, but I also switched the import statements in our main.py. If I then run commit, it not only runs prettier, but it also runs Black and isort. Black was satisfied, isort found that there was an error there and fixed it fast. Prettier fixed it fast as well within the README. The only changes then were in the README, because we already had it nicely sorted before. Hit the second commit. Now the code is nice again. That was pretty much it.

[5:05](https://www.youtube.com/watch?v=08d6zv6zvdM&t=305)
Just a quick shout out to the person who actually brought this tool to us, which was Fabian. It was more than an idea how to always have prettier available in tools without requiring people to install node. He found this tool. Actually it's very nice with having it also run with pipelines and everything. Praise Fabian. With that, I'm open for any questions if there are.

[5:38](https://www.youtube.com/watch?v=08d6zv6zvdM&t=338)
(host) Thank you very much. If you have any questions, then you can either write them in the chat, or you can just ask them straight away. Everyone should be now able to unmute themselves. Do we have any questions? It doesn't seem so. Maxime. Maxime says hello or has a question?

(question) Yes, I have a question. You said that it's in tools. Is it already released, or is it in the coming release.

(answer) In tools, it's already in the release, because actually we run Prettier with pre-commit whenever we dump YAML files or JSON files in tools. For example, when we create the tests for modules, these files are now prettified with Prettier, because before we had the problem that our function actually dumped code, Prettier didn't like. Now we run pre-commit to run Prettier on it. Also the repository itslef has Black and isort. There, it's already in. But for the pipelines, the next release will have it in a template. Then all the pipelines also get it. You just need to run `pre-commit install` to activate it. Modules since yesterday got the got the config in. If you now are in module. Pull. Then run `pre-commit install`, all your modules will automatically run Prettier. Or all changes on modules will automatically be run through Prettier before you can commit. I recommend to do that. If you write modules or subworkflows.

[7:31](https://www.youtube.com/watch?v=08d6zv6zvdM&t=451)
(host) Thank you. Are there any more questions? If not, then I would like to thank our speaker and also the Chan Zuckerberg Initiative for funding the talks. As usual, you can ask more questions if you have any in Slack. This will be uploaded to YouTube. Thank you very much.

(speaker) Bye, everybody.

</details>
