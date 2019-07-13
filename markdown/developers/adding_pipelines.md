---
title: Adding a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
---

# Before you start

So, you want to add a new pipeline to nf-core - brilliant!
Before you start typing, check that you're happy with the following points:

* You're familiar with nf-core and nextflow (see our [introduction docs](/usage/introduction)).
* You're used to working with `git` and [GitHub](https://github.com)
    (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
* The workflow you're thinking of meets the [nf-core guidelines](https://nf-co.re/developers/guidelines).

The main steps involved in adding a new nf-core pipeline covered below are:

1. [Joining the community](#join-the-community)
2. [Creating a pipeline](#create-a-pipeline-from-the-template)
3. [Adding test data](#add-some-test-data)
4. [Adding to the nf-core organisation](#adding-your-pipeline-to-the-nf-core-organisation)
5. [Making your first release](#making-the-first-release)
6. [Updates and new releases](#subsequent-releases)

# Join the community

At its heart, nf-core is a community - to add a pipeline you need to be part of that community!
Please request to join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3))
and introduce yourself on [Slack](https://nf-core-invite.herokuapp.com/) or the
[mailing list](https://groups.google.com/forum/#!forum/nf-core).

It's good to introduce your idea early on so that it can be discussed before you spend lots of time coding.

# Create a pipeline from the template

You'll start by making a new pipeline locally and working with it on your own GitHub account.
Only when it's ready do we move ito the nf-core GitHub organisation.

It's _highly_ recommended to use the nf-core template.
The guidelines for nf-core pipelines are pretty strict, but if you start your pipeline by using the
nf-core template (`nf-core create` - see [the docs](/tools#creating-a-new-workflow))
then your life will be much easier.
This tool does lots of things for you: it gives you the correct file structure and boiler plate code
and also sets up the required `git` infrastructure for you to keep your pipeline in sync in the future.

Even if you already have a working pipeline, it may be easier in the long run to use this this template
and copy over your code in the relevant places.

If you really don't want to use the template it should possible to work without it.
Please see the [manual synchronisation](/developers/sync) documentation.

> Note that workflow names should be all lower-case and contain no punctuation.
> This is to allow consistent names between platforms (eg. GitHub + Docker Hub).

## Push to GitHub

Create a repository on GitHub for your new pipeline under your personal account.

Do this by going to the GitHub website and clicking + then _New Repository_.
Make sure _not_ to initialise it with `README` file - you just want an empty repository.

Once created, copy the URL and add this as a remote to your local git repository
and push your code:

```bash
# Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
# Commit any new code changes since you ran the template
git add .
git commit -m "Starting to build my pipeline"
# Push to GitHub
git push
```

## Set up Travis and Docker Hub

The nf-core pipelines use two additional services that link to GitHub: Travis and Dockerhub.
Using these enables automated testing (Travis) and automatic Docker image builds (Docker Hub).

Just like with GitHub, you can run these on your personal fork for testing purposes.
Once you've merged your pipeline in to the nf-core organisation, we will also set them
up there, but that happens later.

To enable Travis, go to [travis-ci.com](https://travis-ci.com) and log in with your
GitHub credentials. Add a new GitHub repository and select your new pipeline.
The repository should already contain the required `.travis.yml` file, so the next
time you push a commit the tests will be automatically triggered.

The process for Docker Hub is similar, though a little more fiddly. Docker has a number
of  different docker websites (docker hub, docker cloud, docker swarm _etc._), but
they seem to be moving towards consolidating on just Docker Hub. So it's best to use that.

1. Go to [hub.docker.com](https://hub.docker.com) and create an account
2. Create a new repository for your workflow
3. Set your repository to be automatically built from a GitHub repository and link it to your workflow
4. Configure the repo to automatically build whenever you push a new commit to your GitHub repo

Whilst developing your pipeline on your local fork you will need to create automated builds for two docker images
with source set to `master` - one with the `dev` tag and the other with the `latest` tag.
The former will be required for Travis and the latter will be pulled when executing the pipeline locally.

Note: The template assumes that your Docker image is hosted on the nf-core Docker Hub organisation.
To make the pipeline work with your testing image, switch out `nfcore/<PIPELINE_NAME>` for your address
(`username/<PIPELINE_NAME>`). You'll find this in the `container` variable in `nextflow.config` and
`docker` commands in `.travis.yml`.

These will need to be changed back to the defaults before you fork the pipeline to `nf-core`.

## Work on your pipeline

Ok, now you're all set with your own personal nf-core pipeline!
You can now start writing code for real.
Remember to keep running the `nf-core lint` command (see [docs](https://nf-co.re/tools#linting-a-workflow))
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

```bash
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

```bash
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
main nf-core fork. You'll need to ask one of the core nf-core team to help you with this.

## Repository setup

Remember to configure the repository on the GitHub website with the following:

* A description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords!
* Issues enabled, disable Wiki and Projects
* A protected `master` branch that requires review and passing tests
* Write permissions for `nf-core/all` and admin permissions for `nf-core/admin`

## Differences to your own fork

The main difference when working with the main nf-core fork of your workflow is
that tests for pull-requests against the `master` branch will fail. This is because
the `master` branch should only ever contain code from the last release.
Instead, use the `dev` branch for new work and always make pull-requests against
that. Then the tests should pass.

# Making the first release

When the code is stable and ready for a release, make a pull-request from the
`dev` branch to `master` on the nf-core fork. This is a special case and the tests should pass.
Once they do, merge the PR yourself and let the nf-core team know that you're ready.

## Version numbers

When developing the pipeline, the version numbers should be numeric with `dev` at the end.
Use the `nf-core bump-version` command to do this - there are quite a few locations in the
code that need updating and this ensures that it happens in the correct places.
Note that when developing the `:dev` tag should be used for docker containers.

When making a release, version numbers should all be numeric. Use `nf-core lint --release`
when ready - this will check that everything looks correct. You are welcome to use any numeric version number, recommendations are to use [Semantic Versioning](https://semver.org/) as it proved to be a good approach.

## Core pipeline review

Ok - now the tough bit - does your workflow stand up to the scrutiny of the nf-core
team?! Not to worry, we're a friendly bunch. Let us know about the new pipeline,
when you're ready we will create a fake pull-request against the first commit in
the pipeline. This gives the PR review interface showing all code that you've
written. We will go through everything and request and changes that we think are
necessary until you're good to go.

Common things that are flagged at this point are:

* A clear, short but descriptive readme
* Good documentation, especially describing the output files and all parameters
* Pipeline code

We typically tend to have two reviewers for most of the crucial code changes, e.g. adding new major features to an existing pipeline or making an entirely new pipelin release. You can also ping people from the nf-core core team to review your pipelin code by `@`ing them.

## Tagging the release

Once the pseudo-PR is approved, we'll close it and you can [create a new release on GitHub](https://help.github.com/en/articles/creating-releases). Put in a basic changelog entry describing the general functionality at release. You may for example copy the content of your `CHANGELOG` file for that purpose into the GitHub release description. Please use a numeric only release name, to make sure that all pipeline releases follow the same pattern.

The nf-core website and helper tools will automatically detect new releases and update accordingly.

That's it, you're finished! Congratulations!

## Subsequent releases

Once you've made your first release you can continue to work on your fork and make pull-requests
against the `dev` branch on the nf-core repository. Now that we have a stable `master` branch,
there should be reviews of each PR against `dev` before merging.

When ready to make new releases, make sure that the version number is increased and create a
pull-request against `master`. If tests pass, it can be merged and a new release made.

The `master` branch should always have only the commit from the latest release. This is important
because the commit ID is used to reference whether the pipeline is up to date or not.

## Adding new pipeline features to existing pipelines

We are an open and inclusive community, welcoming any contributions to pipelines already present in nf-core. In many cases, the original developers might either not have experience with some new fancy method or simply doesn't have the time to implement everything themselves - so they might be really happy to see you actively contributing!

Basic rules for such contributions:

* Ask in the [Slack](https://nf-core-invite.herokuapp.com/) channel for the specific pipeline whether there is an open issue on the respective pipeline's issue tracker for the feature you're planning to
* If not, create a new issue there, describing the purpose and ideas you have and wait for someone to comment/discuss
* If everyone is happy or there is some consensus in the community, start implementing the feature in your [fork](https://help.github.com/en/articles/fork-a-repo) of the respective pipeline
* Please do not write to multiple channels in the Slack community, rather collect all of the information in a single GitHub issue, which makes it also much easier to follow up on your proposal

### Adding new dependencies to an existing pipeline

Sometimes, especially when adding new features to a pipeline, the dependencies change as well. In such cases, you might want to have an updated Docker Container available before submitting a pull request, in order to have the TravisCI tests run through when testing your updated code. To achieve that, please follow these steps:

* Add *only* the newly required dependencies to the `environment.yml` in the pipeline code
* List this new dependency as something new in the `CHANGELOG`
* Create a Pull Request including only these two changes against the `dev` branch of the pipeline you're working on

This way, a review process will be very fast and we can merge the changes into the `dev` branch, updating the Docker Image for that pipeline automatically. After ~30 Minutes, the Docker Image for that pipeline is then updated, and you can open your Pull Request containing your actual pipeline code changes.

### Release Checklist

As described above, when your pipeline is ready for a (first) release please follow the steps below to ensure that the pipeline meets the requirements of the nf-core community:

* Bump the `dev` branch containing the changes for the release to a release version (e.g. 1.0.0): `nf-core bump-version 1.0.0`
* Check for pipeline dependencies that are out of date and update these accordingly in the `dev` branch: `nf-core lint .` will tell you which ones are outdated via automated API calls to (bio-) conda
* Update the `CHANGELOG`, listing everything that has been added/fixed in this release
* [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) from `dev` to `master` on GitHub after adjusting all of this - waiting for tests to pass
* Link some reviewers (2 are required for merging to `master`) - asking for a final review on your release
* Once approved by two reviewers, merge your PR into `master`
* Go to GitHub and [create a new release for your pipeline](https://help.github.com/en/articles/creating-releases)
  * Please make sure to use strictly numeric release numbers, most people follow [Semantic Versioning](https://semver.org/), e.g. 1.0.0, 1.0.
* Optional: Use a [nice code name](http://www.codenamegenerator.com/) for your pipeline release
* Create your release - tests will automatically run, DockerHub will generate a tagged container for that release.
