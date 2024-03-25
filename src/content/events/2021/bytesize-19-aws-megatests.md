---
title: 'Bytesize 19: Setting up AWS megatests'
subtitle: Gisela Gabernet - QBiC Tübingen, Germany
type: talk
startDate: '2021-09-14'
startTime: '13:00+02:00'
endDate: '2021-09-14'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/2-ekrRsYS00
locations:
  - name: Online
    links:
      - https://youtu.be/2-ekrRsYS00
      - https://doi.org/10.6084/m9.figshare.16621594.v1
      - https://www.bilibili.com/video/BV1w44y1b7bK/
---

This week, Gisela Gabernet ([@ggabernet](http://github.com/ggabernet/)) will tell us all about the nf-core _AWS megatests_‚ and how to set them up in your nf-core pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[1:27](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=87) I’ll be presenting how to set up AWS metatests during this session. These are real-size automated nf-core pipeline tests.

[1:39](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=99) Let’s start by recapitulating how continuous integration of the nf-core pipelines work. For continuous integration, we use GitHub actions, which is a service that is provided by GitHub. This consists of a series of worker instances in which one can run continuous integration tasks. For nf-core pipelines, we have several GitHub actions workflows that define which tasks should be run for continuous integration, and these are listed on the slide.

[2:22](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=142) Just to mention them briefly; one of them is `branch.yml` that takes care of protecting the master branch in the nf-core pipeline repository. It ensures that pull requests coming from external repositories are made to the dev branch and not to the master branch. We also have a GitHub actions workflow that runs pipeline tests, it makes sure that the different test profiles that are added as part of the pipeline run through before integration of the pull request. We also have a couple of linting workflows, some that do code linting and others that lint for trailing white spaces etc. But today, I’m going to focus on two special Github action workflows: The `awsfulltest.yml` and `awstest.yml`, which run the full-size test on AWS batch and the small test data on AWS respectively.

[4:12](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=252) So why did we add those AWS batch tests to the nf-core pipelines? A primary reason was to test that the nf-core pipelines could actually be ported to AWS to be run there without encountering any issues. But we also wanted to test all the nf-core pipelines with full-size data. So the Github action runners that I mentioned before, are limited in their cpu memory and time in the free version of AWS. This wouldn’t allow one to run pipelines with full-size data. So to run pipelines with full-size data, one needs to use something like AWS batch that would allow running those tests. Some additional advantages of having the AWS batch test set up is that we can display results of running full-size data for each of the pipelines on the nf-core website. We can also compare and check that the results remain stable across different releases.

[5:34](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=334) So how do the AWS batch tests actually run? We originally set up the AWS infrastructure that is required to run those tests using [CloudFormation templates](https://docs.opendata.aws/genomics-workflows/orchestration/nextflow/nextflow-overview.html). We needed to adapt them a little, and those were then used to create the AWS infrastructure required for running the batch queues, compute environments etc.

[6:28](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=389) However, recently Seqera has launched Tower Forge, which automatically sets up all of the infrastructure that is needed to run Nextflow pipelines on AWS batch. So we have now migrated and are using Tower Forge to set up all the infrastructure needed to run the pipeline tests. So what does the workflow look like at the moment? We have a Github actions trigger that runs the Github actions workflow that I will show you later. The workflow submits the pipeline job using Tower launch to AWS Batch. The pipeline runs on AWS Batch, pulls the necessary data from the S3 bucket and also publishes the results from the S3 bucket that we have set up. While this is running, we can monitor the progress of the runs via Nextflow Tower. At the end when the pipeline run is finished, we retrieve the results from the S3 bucket and display them on the nf-core website.

[7:53](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=473) So how is this defined in the Github actions workflow? Here is an example of the nf-core AWS full-size test GitHub actions workflow for the nf-core/viralrecon pipeline.

[8:08](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=488) You can see here that it says that the workflow should be triggered on release. But it is also possible to trigger it via a workflow_dispatch button, which I will show you later. What this whole thing does is that it submits the pipeline run using the Tower API. Originally, we were calling the Tower API here, but that was cumbersome, so Phil wrote `tower_action` that allows us to call Tower API via a Tower action and provide the input parameters here in a `.yml` format. This action can also be used outside nf-core pipelines, so feel free to try it out.

[9:24](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=564) So here you can also see all the parameters that are needed to run the jobs and we have most of them stored at secrets as a part of GitHub. We have the pipeline parameters themselves as well as the profiles that are used to run the pipeline. So you don’t really need to modify any of this Github action workflow. This is already set up for you when you use the nf-core template.

[10:01](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=601) The only thing you need to adapt is the test full profile. This is the profile that will be run on AWS when specifying the full test profiles, and it needs to be adapted with all the parameters - including the input sample sheet - to run the full-size test data. For the input sample sheet, this can be set as a part of the pipeline repository. It can also be directly staged on our nf-core AWS megatest S3 bucket and the path provided here.

[10:39](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=639) So let’s summarise what you need to do to set up AWS full-size stats for your pipeline. First, you need to find the full-size data that is suitable for your pipeline. If you don’t have any yet, you can visit the SRA in our repositories and search for some there. Phil has written a tool called [SRA-explorer](https://sra-explorer.info) that allows you to search for data using keywords that might be useful in finding this data. Second, you need to get the data. Check out the [nf-core/fetchngs](https://nf-co.re/fetchngs) pipeline that allows you to download data from SRA or ENA by providing the identifiers. You can also ask a core team member to stage your data on our AWS bucket. Feel free to contact me on Slack if you have any questions. Third, you can add the paths to the data to the sample sheet. Fourth, that’s it! AWS tests will be automatically triggered on release.

[12:56](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=780) Now for a short demo on how that would work. I’ll do my demo on the nf-core/viralrecon pipeline. As I mentioned before, AWS tests are automatically triggered, but in case you'd like to test it prior to release, click on the Actions tab within the pipeline repository and have a look at the actions workflows on the left-hand side of your screen. Running pipelines on AWS is still a bit expensive, so I’m just going to try out a small test for now. So [here](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=829) you see a run workflow button because this workflow has a `workflow_dispatch` event trigger. This can be triggered on any of the branches that are in the repository. So let us select the `dev` branch here and start running the workflow. Now this takes a few seconds to start. But what it does is that it triggers this `github actions` workflow that submits the job to tower. This takes some time.

[14:43](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=883) But I just ran one yesterday, so I can show you how that’s running. So here’s where we monitor the run. We go to Nextflow tower to the AWS megatests workspace within our nf-core organisation. In case you’re a pipeline developer and want to follow your pipeline runs, please contact a member of the core team and we will add you to this workspace.

[15:19](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=917) So let’s have a look and see what it would look like when it’s finished. There are these criteria that highlight the status of the running jobs etc.

[15:36](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=936) If we have a look at the run that I triggered yesterday, it appears that everything ran as it should, and we can also take a look at the allocation of resources for each component of the pipeline.

[15:49](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=949) So once the automated runs are completed, we also display them on the nf-core website. Each pipeline website has a Results tab that holds all the pipeline results; you can take a look at the MultiQC reports for instance.

[16:16](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=976) We request the community to bear in mind that one should work with a reasonable number of samples here because running things on AWS costs us a fair amount of money. We currently have a grant from AWS that covers a limited amount of the costs. You can get an approximate idea of how much the costs were in 2020.

[17:04](https://youtu.be/2-ekrRsYS00?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1024) So that’s it. Please get in touch if you have any questions.

</details>
