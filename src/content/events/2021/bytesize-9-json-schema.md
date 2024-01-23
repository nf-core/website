---
title: 'Bytesize 9: JSON-schema: What, why and how'
subtitle: Matthias Hörtenhuber - Karolinska Institutet, Sweden
type: talk
startDate: '2021-04-20'
startTime: '13:00+02:00'
endDate: '2021-04-20'
endTime: '13:30+02:00'
youtube_embed: https://youtu.be/PU6vAj_7WRM
locationURL:
  - https://youtu.be/PU6vAj_7WRM
  - https://www.bilibili.com/video/BV18h411S76E
  - https://doi.org/10.6084/m9.figshare.14453214.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation.
Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 8: Running pipelines offline

This week, Matthias Hörtenhuber ([@mashehu](http://github.com/mashehu/)) will present: _**JSON-schema: What, why and how**_
This will cover:

- What is a JSON-schema
- How do we use it in nf-core
- How can you create one using `nf-core/tools`

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/PU6vAj_7WRM>
- Slides: <https://doi.org/10.6084/m9.figshare.14453214.v1>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:36](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=36) Welcome everybody, today's topic is JSON schema, what it is, why we chose to use it, and how it can be utilised.

[0:53](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=53) First we need to talk about parameters. In Nextflow you either have them on the command line with `--` in front of it, for example in `main.nf` and in `params.foo`.

[1:12](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=72) When we look through all the pipelines in nf-core, from what I could parse, I found that we have a large number of different parameters, 47 on average but it can go up to 470.

To make this usable, it should be documented so you can read up both as a user and as a developer about what these parameters do.

In the best case, we have a validation and a background so it already throws an error before the tool that uses the parameter. Therefore it is easier to debug.

But as you can imagine, with so many parameters it’s not so easy anymore to manage.

[2:19](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=134) So we therefore looked at JSON schema and thought that this is exactly what we want.

It is a standardised way to describe a data object, and it comes with its own validation library, which takes a load off of our shoulders.

[2:32](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=154) Last summer, we added the `nextflow_schema.json` file to the template, which is a description of all the pipeline parameters in a JSON schema format.

[2:47](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=167) Here you see an excerpt of it, and this is the only thing you should actually see of the JSON because we are of the opinion that you should never have to need to interact with these complex files.

We have tools built around them, so you don’t ever have to touch them.

[3:07](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=189) So how can you make this file?

Or how can you populate it or edit it? It’s simply `nf-core schema build` and then `.` if you’re in your repository.

This then starts the normal wizard, which checks the current file because with a template each pipeline already has some base nf-core parameters in it, and then asks if you want to use the web-builder.

After you click on the blue finish button, it will be sent back here and the parameters will be updated in the nextflow schema json file.

[3:58](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=238) How does it actually look on the website if you run this?

So here you have the parameter schema interface where you can add a parameter.

You scroll there; give it a name, which should be lowercase, you can then also give it an icon, which we render in on all the web-views.

You can give it a description, which doesn’t need to be long, and should in fact be concise enough to fit in the command line help.

If you have more complex things to talk about, we have the help text section, which I will show you shortly.

[4:49](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=289) In the help text, you can write in normal markdown, and it’s interface is something you’ll be familiar with from GitHub.

It also shows you the rendered preview of this markdown.

[5:01](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=301) There are four types of parameters: string, number, integer and a boolean.

For these we have made some special features, but for now say we want a string, give it a default value and decide if it’s a required parameter, which means it always has to be filled out.

We can also say that it should be a hidden parameter, which is used for things you don’t want to be visible for every user.

But it can still be good to interact with it.

[5:41](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=341) If you click on the cog-wheel at the on the side, you see two other features for this parameter.

You can choose for it to be enumerated values and give it different options, or you can give it a pattern, which it will be validated against.

So it is just regular expressions you can use there, and with that you can check for example if an email is valid or if a date or duration is valid.

You can also delete the parameter here if you’ve made a mistake or don’t actually need it.

[6:22](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=380) If we choose that it’s a number, we can get two additional things.

Besides the enumerated values, we can also decide to give it a minimum and a maximum value.

If we are happy there, we can save the changes.

[6:39](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=399) But we at nf-core don’t actually like dangling parameters so we should add it to a group.

[6:46](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=406) A group can have a title, it has a favicon, it can have a description, and also a help text.

You can also hide it and then everything is hidden.

[7:09](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=429) If you now want to quickly put the newly created parameter in, just hit that box or use drag and drop.

When you’re done, you hit finish.

[7:19](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=439) Everything is done here, and on the command line you will see that it stops and wrote the parameters to the new file.

So now we have the updated nextflow json schema.

[7:31](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=451) Why are we using it and where are we using it?

[7:35](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=455) So I stumbled into this because I didn’t want to write documentation in several places for the same thing, because in the main.nf we had the `--help` text where we could define the parameters.

This is now done with the nextflow schema json.

[7:58](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=478) So when you run help, you get this nice output with the groupings intact, the parameter names, parameter types, the description (that’s why you should keep it short), and also the default values that are printed in here.

[8:14](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=494) It’s here, it’s nice, but we also use the same file for the website to render this description of your parameters.

It’s where we also include the icons and the hub text is rendered.

The markdown of the help text is also rendered here, and the default values are shown.

This is how we use it for documentation.

[8:50](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=530) We also use it for parameter validation.

So if you run our pipeline now with our new parameter, the answer and give it 42, it will throw two errors.

The first is that it’s a different type than was expected since it expects a number (since we chose that earlier), instead of a string.

All nf-core pipelines require a parameter, but that was missing here and so that’s the second error.

This is a nice way to quickly check that all the parameters are correct before starting up the whole pipeline.

We also have a launch interface where you just write `nf-core launch <pipeline>`, so either the name of the pipeline or a dot if you’re in the pipeline.

[9:55](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=595) Similar to the how it was for the build, it validates the schema to ensure that everything looks good, and that the parameters match the schema.

Then you are asked if you want to use a web-based or a command line wizard where you can choose the parameters.

I’ll quickly show you later what the web-based interface looks like.

Once we are done with the web-based interface, we are sent back here again. We can decide whether to use the parameters we chose in the browser for this exact run, and run the workflow with our clearly defined parameters.

[10:43](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=643) So what does the web interface look like? Here I have chosen a random pipeline, and you’ve likely seen this `launch` button here, which leads you to a large interface.

You can actually use a different version, and for example I’ve now chosen the development version.

I can launch the development version of the pipeline and the development version of the parameters as well.

[11:14](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=674) Here you see it that it appears to be nicely rendered with a description, with a help text, and also here behind the question mark rendering even complex markdown like tables nicely.

All the parameters here are validated.

[11:30](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=690) Again, I forgot the input parameter, if I add it, I can launch it.

[11:37](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=697) Then we have three options; we can use this ID that takes the parameters we chose for this run, or we can send everything to [Nextflow Tower](https://tower.nf) to use that as an interface, or we can just copy and paste this one into `nf-params.json`.

The latter is actually created if the command line interface is used.

[12:08](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=727) That’s it for how and where we use the nf-core schema, the Nextflow schema json file.

If you want to test your renderings of the markdowns, you can use the website Docker image, and for more documentation, you have the link here in the tools section of the documentation.

If you’d like to discuss how this was implemented or have other questions, join us on the `#json-schema` Slack channel. We are happy to help.

[14:13](https://youtu.be/PU6vAj_7WRM?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=853) You can nest parameters inside.

So according to the definition of json schema, you can have multiple layers of groups inside.

We of course only want one for practical reasons, but we are happy to change this in case the community says they think this would be useful.

Thank you for listening.

</details>
