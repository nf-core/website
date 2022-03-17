---
title: Testing with Gitpod
subtitle: Testing, Code review and Website dev.
---

# Introduction

Gitpod is an open-source developer platform that can be spun up quickly from a Git repository, containing all the programs needed to run your pipelines. It has many purposes, including: testing pipelines, learning Nextflow, reviewing/editing code and live views of edited markdown/HTML for webpage development. In addition, within Gitpod, any change can be pushed to branches for review. 

Gitpod is free (with a GitHub/GitLab/Bitbucket account) for use up to 50hours/month (currently) and based in your browser (reconmmended with Chrome, Firefox or Edge). You can find their extensive documentation [here](https://gitpod.io/). 

All nf-core repos (including the webpage) should have their own working gitpod instances, with all the code required active in each environment.

## Running Gitpod on an nf-core repo

First install the browser extension, follow the instructions [here](https://www.gitpod.io/docs/browser-extension). This adds the green Gitpod button to each Git repository for easy access to the environment. This button simply appends `https://gitpod.io/#` onto the git repo URL (e.g. https://gitpod.io/#https://github.com/nf-core/rnaseq).

Once you click on the green button, you will asked to sign in with either Github, GitLab or Bitbucket. Once you have signed in, you should see the following:

![PNG](/assets/markdown_assets/developers/gitpod/nf-core-gitpod.png)


1. **The terminal** allows you to run all the programs in the repository, for example Nextflow, nf-core tools and docker are installed in the nf-core RNAseq repository.

2. **The sidebar** allows you to customise your environment and perform basic tasks (Copy/Paste, Open files, search, etc.)

3. **The main window** allows you to view and edit files.

Now you can run Nextflow with the chosen nf-core repository. 

For example, for nf-core/rnaseq pipeline we can type the following into the terminal:

    ```console
		nextflow run nf-core/rnaseq \
		-profile test,docker \
		--outdir my_result
    ```

This should run the test data through nf-core rnaseq, using docker with your results in the folder: "my_result". This will take some time to complete.


See the following guides for Gitpod use in different scenarios:

1. [Configuration of a Gitpod environment](gitpod/config.md)
2. [Using git within Gitpod](gitpod/git_in_gitpod.md)
3. [How to test a module with pytest in Gitpod](gitpod/pytest.md)
4. [How to develop Markdown for the nf-core website](gitpod/webdev.md)
5. [Other Gitpod tutorials](gitpod/other.md)



