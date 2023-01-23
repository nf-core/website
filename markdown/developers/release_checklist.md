---
title: Release checklist
subtitle: A step-by-step guide for releasing a nf-core pipeline
---

## Before you release

1. All continuous-integration tests on the `dev` branch should be passing
   - Check the output of `nf-core lint` for warnings and address all that you can
   - Update any software dependencies that are out of date.
     The linting will warn about available updates via automated API calls to (bio-)conda
2. Check that there are no outstanding issues that need to be addressed, especially bug reports.
3. Finalize the description of the pipeline in the GitHub repository and ensure that you remove any "under development" labels in the description. This wording will be used when a new pipeline gets announced on Twitter.
4. If there is a [release milestone](https://help.github.com/en/github/managing-your-work-on-github/about-milestones),
   have a look and see if all issues are closed, or can be resolved
   - It's fine to decide that some things should be postponed until the next release - just move them on to the next milestone

## Steps to release

1. On your fork, bump the version number of the `dev` branch to a release version
   - For example, `1.0.0dev` becomes `1.0.0`
   - Use the `nf-core bump-version` command to make the changes, eg: navigate to the pipeline directory and run `nf-core bump-version 1.0.0`
   - Please make sure to use strictly numeric release numbers
   - Use [Semantic Versioning](https://semver.org/)
2. Run `nf-core lint --release` and check that there are no test failures for release.
3. Check that `CHANGELOG.md` includes everything that has been added/fixed in this release, update the version number above the changes, and optionally add a human-readable release name (e.g. using a [code name generator](http://www.codenamegenerator.com/))
   - We recommend you also add the github tag of the main contributors against each CHANGELOG entry (author, and significant reviewers). This will mean each release on GitHub will display each contributors icons for extra visibility and recognition.
4. [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
5. Once merged, open another PR from the nf-core `dev` branch to the nf-core `master`
   - Make sure that all of the CI tests are passing - this is a special case PR and the tests are different
   - Request PR reviews from at least two people
   - Once approved by two reviewers, merge your PR into `master`
   - And finally delete any label of types: "under development", "under construction" or variants of these
6. Go to GitHub and [create a new release for your pipeline](https://help.github.com/en/articles/creating-releases)
   - NB: Use _exactly_ the same version as in the code (e.g. `1.0.0`) - **do not prefix with v** (e.g. not `v1.0.0`).
   - Optional: Also include your [nice code name](http://www.codenamegenerator.com/) in your pipeline release title (see above with `CHANGELOG.md`)
7. Celebrate! But not too much - you still have a few things left to do...

### Automated events

A number of events are automatically triggered after the pipeline is released:

1. The [nf-core website](https://nf-co.re/pipelines) will be updated automatically with the release information.
2. The [`@nf_core` twitter account](https://twitter.com/nf_core) will send out an automated tweet about the pipeline release within minutes.
3. A [Zenodo DOI](https://zenodo.org/) is automatically generated that provides a persistent means with which to cite the pipeline.

## After release

The last step is to bump up the pipeline version number in the development branch:

1. Bump the version number again on the `dev` branch of **your fork** to a new `dev` version
   - For example, `1.0.0` becomes `1.1.0dev`
   - Use the `nf-core bump-version` command to make the changes, eg: navigate to the pipeline directory and run `nf-core bump-version 1.1.0dev`
2. Update the `CHANGELOG.md` to include a new section for this new version
3. [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.

### After first ever release

After the first release of the pipeline you will need to add the DOI manually into the main `README.md` for the pipeline:

1. Search for your pipeline on Zenodo and find the DOI that allows you to _"Cite all versions"_ of the pipeline.
2. Uncomment the Zenodo-related `TODO` statement in the `Citation` section of the main `README.md` and insert the Zenodo DOI. You should just be able to edit and commit the changes on the `master` branch directly.
3. Add in a badge for the Zenodo DOI at the top of the main `README.md` e.g. [nf-core/atacseq](https://github.com/nf-core/atacseq/blob/fa1e3f8993cd20e249b9df09d29c5498eff311d2/README.md). As with the point above, you should just be able to edit and commit the changes on the `master` branch directly.
