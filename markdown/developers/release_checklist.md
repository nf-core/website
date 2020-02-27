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
  * Use the `nf-core bump-version` command to make the changes, eg: `nf-core bump-version 1.0.0`
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
* Create your release. The tests will run automatically and Docker Hub will generate a tagged container for that release.
  * The nf-core website will automatically update and should automatically send a tweet about the pipeline release within minutes.
  * After about 30 Minutes, the Docker Hub Image for the pipeline will be ready and the released version can be used.

## After release

* Bump the version number again on the `dev` branch to a new `dev` version
  * For example, `1.0.0` becomes `1.1.0dev`
  * Use the `nf-core bump-version` command to make the changes, eg: `nf-core bump-version 1.1.0dev`
* Update the `CHANGELOG.md` to include a new section for this new version
* [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
