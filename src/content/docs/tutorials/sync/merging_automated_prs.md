---
title: Merging automated PRs
subtitle: How to merge automated PRs after a new nf-core/tools release
---

When a new release of tools is created, each pipeline will get an automated pull-request (PR) opened to merge the changes stored in the template into the pipeline.

If there are no merge conflicts on the PR, then that's great!
If you are happy with the changes, feel free to just merge it into the `dev` branch directly.

However, it is quite likely that the PR is quite big with a lot of merge conflicts.
You're going to have to resolve and merge these manually.
Sorry about this, but there's no way around it...

You can either work on the branch created for the template sync to fix the merge conflicts (i.e., on the GitHub web interface), or pull the updates to `TEMPLATE` to your own branch.

Working on your fork is recommended if the merge is not trivial (please make a comment on the automated PR to say that you are working on it though). In this case see the section [Resolving major conflicts](#resolving-major-conflicts) for guidance.

## Resolving minor conflicts

This is the easier route for syncing the template. You can just go to the Pull Requests tab of your repository and open the PR typically named 'Important! Template update for nf-core/tools v1.13.2', which will come from a branch named `nf-core-template-merge-<version>`. This is a modifiable copy of `TEMPLATE`.

At the bottom of the page, resolve the conflicts as guided by GitHub. This should commit to the branch above, and once tests pass you can request reviews from the nf-core community as normal.

## Resolving major conflicts

In the case that there are large conflicts which are unresolvable by the GitHub interface, it is safer and easier to fix these locally in your normal text editor and test on your machine before committing the changes.

The steps we need to do are:

1. Pull the `nf-core/<pipeline>` `TEMPLATE` changes to your fork
2. Resolve the merge conflicts
3. Push these updates to your fork on GitHub
4. Make a PR from your fork to the main nf-core repo

Once you have committed and pushed the updates to your fork and merged these in to the nf-core repository, the automated PR will close itself and show as merged.
You will not need to touch it.

### Pull the changes to your fork

On the command line, go to the directory where you have checked out your fork of the pipeline repository.
Add the nf-core fork as a _git remote_ called `upstream`:

```bash
git remote add upstream https://github.com/nf-core/<pipeline>.git
```

Next, check out a new branch to make these changes in:

```bash
git checkout -b merging-template-updates
```

Finally, pull the `TEMPLATE` branch from the `upstream` repo:

```bash
git pull upstream TEMPLATE
```

### Resolving merge conflicts

You will probably get a tonne of log messages telling you about merge conflicts:

```console
$ git pull upstream TEMPLATE

remote: Enumerating objects: 33, done.
remote: Counting objects: 100% (33/33), done.
remote: Compressing objects: 100% (18/18), done.
remote: Total 33 (delta 15), reused 33 (delta 15), pack-reused 0
Unpacking objects: 100% (33/33), done.
From github.com:nf-core/rnaseq
 * branch            TEMPLATE   -> FETCH_HEAD
   55d617e..2d7814a  TEMPLATE   -> upstream/TEMPLATE
Auto-merging nextflow.config
CONFLICT (content): Merge conflict in nextflow.config
Auto-merging main.nf
CONFLICT (content): Merge conflict in main.nf
Auto-merging environment.yml
CONFLICT (content): Merge conflict in environment.yml
...
```

If you look at the current status, you will see the files that have merge conflicts that need resolving _(Unmerged paths)_:

```console
$ git status

On branch merging-template-updates
You have unmerged paths.
  (fix conflicts and run "git commit")
  (use "git merge --abort" to abort the merge)

Changes to be committed:

    modified:   .github/ISSUE_TEMPLATE/bug_report.md
    modified:   .github/ISSUE_TEMPLATE/feature_request.md
    modified:   .github/markdownlint.yml
    modified:   .gitignore
    new file:   bin/markdown_to_html.py
    deleted:    bin/markdown_to_html.r
    deleted:    conf/awsbatch.config

Unmerged paths:
  (use "git add/rm <file>..." as appropriate to mark resolution)

    both modified:   .github/CONTRIBUTING.md
    both modified:   .github/PULL_REQUEST_TEMPLATE.md
    both added:      .github/workflows/branch.yml
    both added:      .github/workflows/ci.yml
    both added:      .github/workflows/linting.yml
    deleted by them: .travis.yml
    both modified:   CHANGELOG.md
    both modified:   CODE_OF_CONDUCT.md
    both modified:   Dockerfile
    both modified:   README.md
    both modified:   assets/multiqc_config.yaml
    both modified:   bin/scrape_software_versions.py
    both modified:   conf/base.config
    both modified:   conf/igenomes.config
    both modified:   conf/test.config
    both modified:   docs/output.md
    both modified:   docs/usage.md
    both modified:   environment.yml
    both modified:   main.nf
    both modified:   nextflow.config
```

You now need to go through each of these files to resolve every merge conflict.
Most code editors have tools to help with this, for example [VSCode](https://code.visualstudio.com/docs/editor/versioncontrol#_merge-conflicts) have built-in support.

Be careful when resolving conflicts.
Most of the time you will want to use the version from the `TEMPLATE` branch,
but be aware that some of this new template code may need to be customised by your pipeline.
In other words, you may need to manually combine the two versions in to one new code block.

If you have any doubts, ask for help on the nf-core Slack.

### Pushing the resolved changes to your fork

When all merge conflicts have been resolved and all files are staged, you can commit and push these changes as with any other new code:

```bash
git commit -m "Merged changes from nf-core template"
git push --set-upstream origin merging-template-updates
```

### Merging to the nf-core repository

Once the changes are on your fork, you can make a pull request to the main nf-core repository for the pipeline.
This should be reviewed and merged as usual.
You should see in the commit history on the PR that there is a commit by the @nf-core-bot user, with the same commit hash found in the automated `TEMPLATE` PR.

Once your fork is merged, the automated PR will also show as merged and will close automatically.
