---
title: 'Bytesize: nf-validation'
subtitle: Júlia Mir Pedrol, QBIC and Nicolas Vannieuwkerke, Center for Medical Genetics Ghent
type: talk
start_date: '2023-06-06'
start_time: '13:00+02:00'
end_date: '2023-06-06'
end_time: '13:30+02:00'
youtube_embed: https://www.youtube.com/watch?v=rr9FTlQayIE
location_url:
  - https://www.youtube.com/watch?v=rr9FTlQayIE
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-validation

This week Júlia Mir Pedrol ([@mirpedrol](https://github.com/mirpedrol)) and Nicolas Vannieuwkerke ([@nvnieuwk](https://github.com/nvnieuwk)) are presenting nf-validation, the soon to be released plugin to validate the input parameters and sample sheets for Nextflow pipelines.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1)
Hello, everyone, and welcome to today's bytesize talk. I'm happy to introduce to you Julia and Nicolas, and they're going to talk about nf-validation. I'm handing over now to you, Julia.

[0:20](https://www.youtube.com/watch?v=rr9FTlQayIE&t=20)
Thank you. Hello, everyone. We're going to explain this new plugin that we implemented in Nextflow. It's called nf-validation, and we use it for pipeline parameter validation, and for this we use JSON schema. First of all, before starting, why is it important to validate parameters? You may know that Nextflow pipelines can accept different parameters, either through command line or through all the config files, and this is not validated by Nextflow. If, for example, your pipeline expects a string and the user provides a number, all of your pipelines will run until this value is used, and then it will fail. That's why it's important to have some previous steps to validate parameters and avoid possible errors. This was already implemented in nf-core using this JSON schema, and actually all the nf-core pipelines have these validation steps in the template, because if, as a pipeline developer, you would have to validate these things manually, it would be a huge chunk of code.

[1:45](https://www.youtube.com/watch?v=rr9FTlQayIE&t=105)
We use JSON schemas, as I said, and this JSON schema looks something like that. Here you describe all the parameters of your pipeline. It has some formatting, and then under these definitions we have groups, because you can organize your parameters by, for example, input parameters and at least organize them in different groups. Then inside properties you have, for example, in this case, foo, which should be string, or bar, which should be string. This file can get very long, the advice is to never edit it by hand. In nf-core, I think there's another bytesize talk about that, but in nf-core you have this command, `nf-core schema build`, which will open a web tooling, which helps edit this JSON file, this JSON schema, and it's like a drag and drop, so it's very easy to edit it and you don't need to be careful with the formatting and so on.

[3:14](https://www.youtube.com/watch?v=rr9FTlQayIE&t=194)
Then another thing that is new from this plugin is that... This JSON schema can be used for different things, for example, we also use it in the nf-core website, but it can also validate other kinds of files, for example, sample sheets, which usually are used in pipelines to provide inputs, so it's usually a CSV or a TSV file, and where you have your sample ID. If you provide files, you can have each column providing one file and maybe some metadata from samples or things like that. You can also have a JSON schema to validate this sample sheet. The format is more or less the same as the one that I already showed, the structure is tiny bit different, but you also have properties, and inside properties you would have the name of each column in your CSV or TSV. It can also validate YAML files. In this case you will have the name of every entry, and then you can also have type, and you can validate different things. For example in the case of being a string, you will validate that the provided value is a string, or you can also provide a pattern if it has to end with .fasta or things like that.

[5:03](https://www.youtube.com/watch?v=rr9FTlQayIE&t=300)
We'll go now, this was a little bit fast, but so I think we have another bytesize about JSON schema which are more in detail, but for the time now, I'm going to talk about the nf-validation, the plugin itself. This plugin takes all the code that was started in nf-core. If you have checked the nf-core template at some point, this is how the pipeline template looks, and you have here a lib directory. Then here we store some group code, and for example this file is the one that validates Nextflow parameters. This was taken from the nf-core template, and based on that we started the development of this plugin. How to use it is very easy. Like all Nextflow pipelines, you can add in your Nextflow .config, these plugins, and then you add the name of the plugin you want to use, and the latest version. With this, that's all what you need, then this will be installed with your Nextflow, and then in this case it contains different functions that can be imported in your main.nf or in your Nextflow script, and you only need to include the name of the function that then you can use in your script, and from plugin nf-validation.

[7:09](https://www.youtube.com/watch?v=rr9FTlQayIE&t=429)
Then these functions we have here, we have different ones. I will quickly go through them as a summary. We have paramsHelp(), which is used to print a help message for a pipeline, so just a show. You would use... now I'm using launch.sh because the latest version is not released, so that's running my local copy, but usually it would be `nextflow run`, and then if you have this in your nextflow.config you don't need the remote, the name of your pipeline, and then we can run help. This uses a JSON schema that I talked about to print the help message of the pipeline. If it's not working... yes, perfect. Here you see the help message with the usual command, and then the parameters, those are the sections where they are organized, and then you see the name and some description, also the type of value. Then we also have paramsSummaryLog() and paramsSummaryMap(). These two work very similar, and they are used to print. Usually when you run a pipeline, in nf-core at least we print a summary of the parameters that change from the default at the beginning of every run, in case a user needs to check what they provided. This is generated with this function, paramsSummaryLog(), which provides this list of parameters in text format, and paramsSummaryMap() works exactly the same but instead of returning a text format, it returns a map.

[10:09](https://www.youtube.com/watch?v=rr9FTlQayIE&t=609)
Then we also have validateParameters(), which is maybe the most important here, which is the one that does the actual validation of the parameters. In your main.nf, you can use the function validateParameters(), and if you use this function before starting the execution of the workflow, it will fail in case there's some error before starting all the execution. For example here it says the parameter that you provided called input, it's sample sheet text and it doesn't match the pattern csv, tsv or yaml, and also it's a file that doesn't exist, it's also validating that this file should exist.

[11:11](https://www.youtube.com/watch?v=rr9FTlQayIE&t=671)
I'm gonna show as an example how this looks, so that's the current template without using the plugin that we have in nf-core, and as you see we use this chunk of code which is initialising and also validating all the parameters, and then here I have the same template but modified in order to use the plugin. Here I imported the functions and instead of... before I had this initialise, which was using all the code inside lib, in this case this has been modified and we don't have any more the nf-core schema.groovy, and I have the code to print a help message and here the function to validate parameters. If I run this pipeline again, the test for example, it should validate all the parameters and now we will see first the summary of parameters that I mentioned before. As you see, because I didn't provide the outdir parameter, which is required, I get this error before starting any execution. That's the description of the parameters that are different from before, for example you can see which input file you provided, and then now the validation passed and our pipeline started.

[13:44](https://www.youtube.com/watch?v=rr9FTlQayIE&t=824)
And then the last function that we have is fromSamplesheet(), which is reading the input sample sheet and creating a channel, and I will leave this for the end because Nicolas worked on that so he will explain about this. Also a new thing that we have now with this plugin: you can have schemas inside schemas. What does this mean? In your original Nextflow schema file, for every parameter which is a file, you can have this new key called schema, and this one references to a path of another JSON schema, in this case it's a JSON schema which will validate the input sample sheet. This will also, now you'll see it when Nicolas explains more in detail, so it's automatically whenever it detects that there's this schema key in a parameter, it will try to read the file provided by this parameter and then validate it using this JSON schema.

[15:12](https://www.youtube.com/watch?v=rr9FTlQayIE&t=912)
I also have a different example here for RNA-seq. If we see the main code, that's exactly the same that I showed before where I import the functions, also bring the help message and validate parameters. RNA-seq was one of the first pipelines that got this input schema, it was just like a proof of concept. We started implementing this some time ago, and now we have it implemented with a plugin so all pipelines can use it. Here you have the columns of your sample sheet, in this case sample ID, FASTQ1, FASTQ2 and strandedness. This will automatically validate the content of the input. I was gonna try to run this pipeline but maybe I'm talking too much and it's a bit long so that's it, that you know that you can now automatically validate, and that's what works for sample sheets but also for any other TSV, CSV or YAML file that you would like to validate. It doesn't have to be the input specifically, it will validate any of these files. So now I will hand over to Nicholas if he wants to show this fromSamplesheet().

[17:03](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1023)
Yes, let me share my screen. Okay so I'm gonna show you a real quick example of how to use the fromSamplesheet() function. As you can see, I have a simple pipeline here which validates parameters and then converts the input parameter sample sheets to a channel. I'm going to run it real quick, as you can see, I also use the launch.sh bash script because the version is not released yet. If I'm running this, you can see I get some outputs which is the channel inputs. This output has been made from the sample sheet CSV as you can see, it has the name, surname, the likes and pictures from certain persons, for example the first line Harry Potter, the full path to a text file which has his likes in it and a full path to a directory which has pictures in it which correspond to his likes. I use this sample sheet to validate the schema, to validate the sample sheets.

[18:15](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1095)
As you can see, all properties are inside of an items section. You can see the name, surname, a hidden ID number which isn't in the sample sheets but you see if it isn't in the sample sheet it will automatically go to null. Then it has the likes which has a format file path and it checks if the file exists and also pictures which is directory which also contains the key dependent required. If likes is not given but pictures is given, the sample sheet validation will fail. I'll show an example of this. For example, if I remove Harry Potter's likes and rerun the codes, you'll see an error which says that the likes fields should be defined when pictures is specified. Then it also shows which fields are not defined because you can also add way more fields to it and I think it's a very nice error message. You can also specify the unique which will take a boolean or a list. If it's boolean it will only look at the field itself, so all names should be unique if it's true. If you give a list with, for example, surname, all fields should be unique together with the surname. For example, I can't specify Harry Potter twice so you can see this. It also gives the error for the likes fields which is not specified. But as you can see the combination of name with field surname needs to be unique and you see which combination is the one that clashes with it not being unique. Okay, so this is a real small example of how the fromSamplesheets() works.

[20:15](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1215)
One small thing to note is that the unique and dependent required field parameters actually only validates if you run the fromSamplesheet()s, because these are specific for the sample sheet conversion and won't be validated using validateParameters(). All the other schema fields will be validated using validatePrameters(). One other nice thing with fromSamplesheet() is that it will create meta fields which are immutable from the start. I have a bit of code here to show you this. If I try to change the name of every character to Voldemort it will fail because it cannot change the value in a meta field. Of course I have to make sure my sample sheet passes first. As you can see, you cannot put items into an immutable map. This can cause problems though in some pipelines which are already built around this concept and so you can disable it using the optional key immutable meta by defining FALSE. The default of this is TRUE as you can see.

[21:34](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1294)
If I define this and run it again I will be able to change the name of every character to Voldemort, or it should do that. Apparently it does not. I don't know why. You can also do it with a parameter let's see if that works. It does not work okay so normally that should work I think I made a typo somewhere or something. You can find this all in the documentation of the validation plugin. You can also specify this fromSamplesheet() function from default. Go to assets/schema_input.json/schema file to convert sample sheets to a channel, which can also specify which schema to use by using schema, then path to schema. It's weird that it's not working. I'll try it again. No. Okay, weird. That's it for the fromSamplesheet()s conversion. Any questions?

[22:59](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1379)
(host) Thank you very much.

(speaker) I just wanted to share the last last slide.

(host) I'm so sorry!

(speaker) No no, it's fine. Just to share the point... what am I sharing? Just a quick thanks to Phil and Kevin. Kevin started this code in nf-core and also to everyone who contributed on nf-core to this either by testing or reviewing documentation, especially. Just to share some important things that you may want to check. The repo of nf-validation is in Nextflow and here you have these documentation, that I have been also using. It's a very nice documentation and pretty extensive. We have this Slack channel which is shared in nf-core and Nextflow called nf-validation. The last thing to mention is that this will be coming soon in the nf-core template in the next release. The parameter validation and also optional, not mandatory but optional, obtaining this input channel with fromSamplesheet(). That's everything, thank you.

[24:29](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1469)
(host) Thank you again. There was someone who had a question I think. You can now unmute yourself and also start the video if you have a question.

(question) Hi, thank you for this presentation. It's very nice. I was wondering, Julia, will you be adding a schema build command for the sample sheet, too?

(speaker) An schema, what? Sorry, can you say that again?

(question cont.) An nf-core schema build command for the sample sheet

(answer) Yes, exactly. that's not existing right now. Now we have this tooling to create the Nextflow schema for parameters but not for sample sheets, but there's the plan to add it and highly probable to move all this tooling out of nf-core and make it also as a standalone. I guess I can't estimate a date, because that's quite a bit of work. We'll see.

(question cont.) All right, thank you.

[25:49](https://www.youtube.com/watch?v=rr9FTlQayIE&t=1549)
(host) Are there any more questions from the audience? It doesn't seem so. In that case, I would like to thank both of you and of course the audience for listening and as usual the Chan Zuckerberg Initiative for funding our bytesize talks. This will be the last bytesize talk before our summer break and we will let you know when we commence after summer, so thank you very much everyone.

</details>
