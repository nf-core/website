---
title: Release checklist
subtitle: A step-by-step guide for releasing a nf-core pipeline
---

## Before you release

* All continuous-integration tests on the `dev` branch should be passing
  * Check the output of `nf-core lint` for warnings and address all that you can
  * Update any software dependencies that are out of date.
    The linting will warn about availble updates via automated API calls to (bio-)conda
* Check that there are no outstanding issues that need to be addressed, especially bug reports
* If there is a [release milestone](https://help.github.com/en/github/managing-your-work-on-github/about-milestones),
  have a look and see if all issues are closed, or can be resolved
  * It's fine to decide that some things should be postponed until the next release - just move them on to the next milestone

## Steps to release

* Bump the version number of the `dev` branch to a release version
  * For example, `1.0.0dev` becomes `1.0.0`
  * Use the `nf-core bump-version` command to make the changes, eg: `nf-core bump-version <path_to_cloned_pipeline> 1.0.0`
  * Please make sure to use strictly numeric release numbers
  * If in doubt, use [Semantic Versioning](https://semver.org/) as a guide
* Run `nf-core lint` with the `--release` flag and check that there are no test failures.
* Check that `CHANGELOG.md` includes everything that has been added/fixed in this release
* [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
* Once merged, open another PR from the nf-core `dev` branch to the nf-core `master`
  * Make sure that all of the CI tests are passing - this is a special case PR and the tests are different
  * Request PR reviews from at least two people
  * Once approved by two reviewers, merge your PR into `master`
* Go to GitHub and [create a new release for your pipeline](https://help.github.com/en/articles/creating-releases)
  * Optional: Use a [nice code name](http://www.codenamegenerator.com/) for your pipeline release
* Create your release.

## After release

A number of events are automatically triggered after the pipeline is released:

* A tagged container for the release will be built in the [nf-core Docker Hub](https://hub.docker.com/orgs/nfcore/repositories) account. This will take ~30-45 minutes to build after the release.
* The [nf-core website](https://nf-co.re/pipelines) will be updated automatically with the release information.
* The [nf-core Twitter](https://twitter.com/nf_core) account will send out an automated tweet about the pipeline release within minutes.
* A [Zenodo DOI](https://zenodo.org/) is automatically generated that provides a persistent means with which to cite the pipeline.

After the first release of the pipeline you will need to add the DOI manually into the main `README.md` for the pipeline:

* Search for your pipeline on Zenodo and find the DOI that allows you to _"Cite all versions"_ of the pipeline.
* Uncomment the Zenodo-related `TODO` statement in the `Citation` section of the main `README.md` and insert the Zenodo DOI. You should just be able to edit and commit the changes on the `master` branch directly.
* Add in a badge for the Zenodo DOI at the top of the main `README.md` e.g. [nf-core/atacseq](https://github.com/nf-core/atacseq/blob/fa1e3f8993cd20e249b9df09d29c5498eff311d2/README.md). As with the point above, you should just be able to edit and commit the changes on the `master` branch directly.

Finally, don't forget to bump up the pipeline version number in the development branch:

* Bump the version number again on the `dev` branch to a new `dev` version
  * For example, `1.0.0` becomes `1.1.0dev`
  * Use the `nf-core bump-version` command to make the changes, eg: `nf-core bump-version 1.1.0dev`
* Update the `CHANGELOG.md` to include a new section for this new version
* [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
