# Before you start
So, you want to add a new pipeline to nf-core - brilliant!
Before you start typing, check that you're happy with the following points:

* You're familiar with nf-core and nextflow (see our [introduction docs](/usage_docs)).
* You're used to working with `git` and [GitHub](https://github.com)
    (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
* The workflow you're thinking of meets the [nf-core guidelines](/guidelines#guidelines-for-nf-core-pipelines).

# Join the community
At its heart, nf-core is a community - to add a pipeline you need to be part of that community!
Please request to join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3))
and introduce yourself on [Slack](https://nf-core-invite.herokuapp.com/) or the
[mailing list](https://groups.google.com/forum/#!forum/nf-core).

It's good to introduce your idea early on so that it can be discussed before you spend lots of time coding.

# Create a pipeline from the template
The guidelines for nf-core pipelines are pretty strict, but if you start your pipeline by using the
nf-core template (`nf-core create` - see [the docs](http://localhost:8888/tools#creating-a-new-workflow))
then your life will be much easier.
This tool does lots of things for you: it gives you the correct file structure and boiler plate code
and also sets up the required `git` infrastructure for you to keep your pipeline in sync in the future.

Even if you already have a working pipeline, it may be easier in the long run to use this this template
and copy over your code in the relevant places.

If you really don't want to use the template it should possible to work without it.
Please see the [manual synchronisation](/sync) documentation.

> Note that workflow names should be all lower-case and contain no punctuation.
> This is to allow consistent names between platforms (eg. GitHub + Docker Hub).

## Push to GitHub
Create a repository on GitHub for your new pipeline under your personal account.
Make sure _not_ to initialise it with `README` file - you just want an empty repository.

Once created, copy the URL and add this as a remote to your local git repository
and push your code:

```bash
# Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
# Commit any new code changes
git add .
git commit -m "Starting to build my pipeline"
# Push to GitHub
git push
```

> NB: We hope to automate these steps with the `nf-core create` command soon.

## Set up Travis and Docker Hub
The nf-core pipelines use two additional services that link to GitHub: Travis and Dockerhub.
You'll want to set both of these up for your new pipeline to enable automated
testing (Travis) and automatic Docker image builds (Docker Hub).

To enable Travis, go to [travis-ci.org](https://travis-ci.org) and log in with your
GitHub credentials. Add a new GitHub repository and select your new pipeline.
The repository should already contain the required `.travis.yml` file, so the next
time you push a commit the tests will be automatically triggered.

The process for Docker Hub is similar, though made a little complicated by the number
of different docker tools (docker hub, docker cloud, docker swarm _etc._).
Docker Hub and Docker Cloud use the same back-end, so it doesn't really matter which
website you use. Docker Cloud has a nicer interface though, so we recommend that.

1. Go to [cloud.docker.com](https://cloud.docker.com) and create an account
2. Create a new repository for your workflow
3. Set your repository to be automatically built from a GitHub repository and link it to your workflow
4. Configure the repo to automatically build whenever you push a new commit to your GitHub repo

Whilst developing your pipeline on your local fork you will need to create automated builds for two docker images
with source set to `master` - one with the `dev` tag and the other with the `latest` tag.
The former will be required for Travis and the latter will be pulled when executing the pipeline locally.

> NB: The default name (`nfcore/<PIPELINE_NAME>`) for the Docker image in the pipeline template
will need to be replaced i.e. `container` variable in `nextflow.config` and `docker` commands in `.travis.yml`.
These will need to be changed back to the defaults before you fork the pipeline to `nf-core`.

## Work on your pipeline
Ok, now you're all set with your own personal nf-core pipeline!
You can now start writing code for real.
Remember to keep running the `nf-core lint` command (see [docs](http://localhost:8888/tools#linting-a-workflow))
to make sure that your workflow passes all of the nf-core tests.
The automated tests on Travis also run this so you should get an email if something breaks.

# Add some test data
Whilst the linting tests are good, they're not sufficient by themselves.
It's also good to get Travis to actually run your pipeline on a minimal dataset.
Currently, we don't usually check the results that are produced, but it often catches
syntax errors and other serious problems that cause nextflow to exit with an error.

## Putting the test data on GitHub
To avoid bloating the workflow, we don't keep test data in the same repository as
nf-core workflows.
Instead, we use the dedicated [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.

To set this up, make a fork of that repo to your personal account.
Clone the repository and check out a new branch for your workflow:

```
git clone https://github.com/YOUR_USERNAME/test-datasets.git
cd test-datasets
git checkout -b MY_WORKFLOW
```

Now add your test data files - note that they must be **very small**.
GitHub has quite a low file size limit, and the Travis jobs will time out with anything
that's not tiny. We typically use PhiX / Yeast / part of a chromosome as a reference
and aggressively subsampled input data.

Once added, push these new files to GitHub:
```bash
git add .
git commit -m "Added test data for MY_WORKFLOW"
git push --set-upstream origin MY_WORKFLOW
```

Finally, make a pull-request against the main [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository with your files.
You want this repo to also use a branch with the name of your workflow, so first go
to the [repository GitHub web page](https://github.com/nf-core/test-datasets/) and create
this new branch using the UI there.
Once created, you can open a pull request and select this as the target branch.

If in doubt, ask for help!
([Slack](https://nf-core-invite.herokuapp.com/) or [mailing list](https://groups.google.com/forum/#!forum/nf-core))

## Setting up a test workflow
Now that your test data is hosted on the web, you can set up a `test` config profile in your
workflow that points to it.
In fact, the `test` profile should already exist if you've used the template.
Switch out the example URLs for the ones you added (view the files on GitHub and click 'Raw' to get the URL).

Add any other required parameters so that running the pipeline runs with as few extra
flags as possible. Note that the `test` profile can be combined with other profiles such as `docker`
or `conda`, so your config should not specify a hardware environment.

Have a go at running the pipeline and see if it works:
```
nextflow run MY_WORKFLOW -profile test,docker
```

Note that if you do need to adjust this `nextflow run` command, you'll need to update it
in the `.travis.yml` config file too.


# Adding your pipeline to the nf-core organisation
Ok, so you're essentially finished. Your pipeline is written, the tests pass and
you're ready to add your workflow to nf-core.

First, fork your workflow repository to the nf-core GitHub organisation by
clicking 'Fork' at the top of the GitHub webpage. If you don't see nf-core
as an option, please ask one of the nf-core administrators to do this for you.

Once forked, the [nf-core website](https://nf-co.re) will automatically update
to list your new pipeline.

## Setting up Travis and Docker Hub
Just as with your own fork, Travis and Docker Hub need to be set up for the
main nf-core fork.
If in doubt, please ask an nf-core administrator to help you with this.

## Repository setup
You or someone else should set up the new repository with the following:
* A description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords!
* Issues, no Wiki or Projects
* A protected `master` branch that requires review and passing tests
* Write permissions for nf-core/all and admin permissions for nf-core/admin

## Core pipeline review
Ok - now the tough bit - does your workflow stand up to the scrutiny of the nf-core
team?! Not to worry, we're a friendly bunch. Let us know about the new pipeline,
when you're ready we will create a fake pull-request against the first commit in
the pipeline. This gives the PR review interface showing all code that you've
written. We will go through everything and request and changes that we think are
necessary until you're good to go.

## Differences to your own fork
The main difference when working with the main nf-core fork of your workflow is
that tests for pull-requests against the `master` branch will fail. This is because
the `master` branch should only ever contain code from the last release.
Instead, use the `dev` branch for new work and always make pull-requests against
that. Then the tests should pass.

# Making a release
When the code is stable and ready for a release, make a pull-request from the
`dev` branch to `master` on the nf-core fork. This is a special case and the tests should pass.
Once they do, merge the PR and create a new release on GitHub.

The nf-core website and helper tools will automatically detect new releases and update accordingly.

That's it, you're finished! Congratulations!
