---
title: Synchronisation
subtitle: Learn how nf-core pipelines are automatically kept up to date with community standards.
---

# Introduction

To keep all the nf-core pipelines up-to-date with the latest version of the community standards, we have implemented a synchronisation tool.
This ensures that updates to the community standards are propagated to all nf-core pipelines.

There are three topics covered in this documentation page:

1. [Merging automated PRs](#merging-automated-prs)
2. [Manual synchronisation](#manual-synchronisation)
3. [Setting up a pipeline for syncing retrospectively](#setting-up-a-pipeline-for-syncing-retrospectively)
4. [Fixing a broken `TEMPLATE` branch](#fixing-a-broken-template-branch)

### How template synchronisation works

The `nf-core` helper tools have a subcommand for synchronising a pipeline with the nf-core template (`nf-core sync`).
Although this can be run manually, it is usually only used by the GitHub Actions automation:
when a new version of `nf-core/tools` is released it runs for all nf-core pipelines and automatically opens pull-requests (PRs) with the necessary changes required to update every pipeline.
These pull requests then need to be manually resolved and merged by the pipeline maintainers.

Behind the scenes, this synchronisation is done by using `git`.
Each repository has a special `TEMPLATE` branch which contains only the "vanilla" code made by the `nf-core create` tool.
The synchronisation tool fetches the essential variables needed to recreate the pipeline and uses this to trigger a `nf-core create --no-git` command with the latest version of the template.
The result from this is then compared against what is stored in the `TEMPLATE` branch and committed. During an automated sync, a copy of the `TEMPLATE` branch called `nf-core-template-merge-<version>` will be made (to avoid `dev` history to end up in `TEMPLATE` branch after solving merge conflicts), and a PR from this new branch will be opened against your `dev`.
When merging from the `nf-core-template-merge-<version>` branch back into the main `dev` branch of the pipeline, `git` should be clever enough to know what has changed since the template was first used, and therefore, it will only present the relevant changes.

For this to work in practice, the `TEMPLATE` branch needs to have a shared `git` history with the `master` branch of the pipeline.
The `nf-core create` command initially does this by enforcing a first commit to the `master` branch before any development has taken place.
If the pipeline _was not_ created by the `nf-core create` command, this has to be set up manually.
For instructions on this, see [Setting up a pipeline for syncing retrospectively](#setting-up-a-pipeline-for-syncing-retrospectively).

# Merging automated PRs

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

# Manual synchronisation

There are rare cases, when the synchronisation needs to be triggered manually,
i.e. it was not executed during an `nf-core/tools` release on Github, or when you want to perform a targeted sync.

Note that automated PR system is only applicable to official nf-core pipelines, homemade pipelines based on nf-core standards/modules created with `nf-core create` have to be updated following this manual synchronisation procedure.

You can do so by running the `nf-core sync` command:

```bash
cd my_pipeline
git checkout dev # or your most up to date branch
nf-core sync -d .
```

Note that the `sync` command assumes that you have a branch called `TEMPLATE`, so you may need to pull this from the upstream nf-core repository if you are working on a fork:

```bash
git remote add upstream https://github.com/nf-core/PIPELINE.git
git checkout --track upstream/TEMPLATE
```

Remember to go back to your `dev` branch as above before running `nf-core sync`.

Much of the merging process should then be the same as described above with the automated pull requests.

In case of manual synchronisation of a homemade pipeline and if you want to have a PR opened to your `dev` branch, you can use this `nf-core sync` template command:

```bash
nf-core sync \
   --dir [pipeline_dir]
   --from-branch dev \
   --pull-request \
   --username [GitHub_username] \
   --github-repository [GitHub_pipeline_URL]
```

# Setting up a pipeline for syncing retrospectively

This section describes how to set up a correct TEMPLATE branch in the case your pipeline was not created with a TEMPLATE branch from the beginning. If you created a pipeline with the `nf-core create` command, you should be all ready to go and can skip this step. Otherwise proceed with caution. It is probably a good idea to make sure you have all your local changes pushed to github and you could even make a local backup clone of your repository before proceeding.

You should also consider the option to restart your pipeline project by running the `nf-core create` command and simply copy in the modifications you need into the newly created pipeline.

### Step-by-step procedure

This walkthrough assumes that you are working directly with the head `nf-core` fork of the pipeline.
It is possible (and potentially safer) to do this on your own fork instead, it's up to you.

First clone your pipeline into a new directory (in case we mess things up):

```bash
mkdir TMPDIR
cd TMPDIR
git clone https://github.com/nf-core/YOURPIPELINENAME.git
```

Then create the new TEMPLATE branch and _delete all your files_ in order to have a completely empty branch:

```bash
cd pipeline_root_dir
git checkout --orphan TEMPLATE && git rm -rf '*'
```

Make sure your branch is completely empty by checking the status of `git status`:

```bash
$ git status
On branch TEMPLATE

No commits yet

nothing to commit (create/copy files and use "git add" to track)
```

Regenerate your pipeline from scratch using the most recent template:

> Make sure you are within your pipeline root directory before running these commands.

```bash
nf-core create --no-git
```

If your pipeline already has versioned releases (eg. you are not currently on `1.0dev`),
then specify the version number that you are currently on:

```bash
nf-core create --no-git --version 1.3dev
```

> The version you choose should match the branch that you intend to merge with.
> If you already have a release, you should probably be merging in to `dev` eventually,
> so use the version number specified there.

Follow the prompts to fill in the pipeline name, description and author(s).
Make sure that you take the exact text that you already have already used in your pipeline's
`nextflow.config` file (`manifest.name` etc.), if these have already been written.

This creates a new directory `YOURPIPELINENAME` with the template pipeline files in it.
Now move these files into your root git directory:

```bash
mv nf-core-YOURPIPELINENAME/* .
mv nf-core-YOURPIPELINENAME/.[!.]* .
rmdir nf-core-YOURPIPELINENAME
```

Now make sure the newly created files are in the correct place. It should look similar to this:

```console
$ git status
On branch TEMPLATE

No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)

  .editorconfig
  .gitattributes
  .github/
  .gitignore
  .markdownlint.yml
  CHANGELOG.md
  CITATIONS.md
  CODE_OF_CONDUCT.md
  LICENSE
  README.md
  assets/
  bin/
  conf/
  docs/
  lib/
  main.nf
  modules.json
  modules/
  nextflow.config
  nextflow_schema.json
  subworkflows/
  workflows/

nothing added to commit but untracked files present (use "git add" to track)
```

If it all looks good, then commit these files:

```bash
git add .
git commit -m "Initial template commit"
```

For the nf-core bot to be able to access your `TEMPLATE` branch, you need to push it to the upstream repository ([https://github.com/nf-core](https://github.com/nf-core)).

```bash
git push --set-upstream origin TEMPLATE
```

### Merge TEMPLATE into main branches

The only remaining step is unfortunately a rather tedious one.
You have to merge the `TEMPLATE` branch into your main pipeline branches, manually resolving all merge conflicts.

If your pipeline is in early development, you can do this with `master` branch directly. If not, it's better to do this
in a branch and then you can make a pull-request to `dev` / `master` when ready.

```bash
git checkout dev
git checkout -b template_merge
git merge TEMPLATE --allow-unrelated-histories
```

You can try extra flags such as `-Xignore-space-at-eol` if you find that the merge command shows entire files as being new.

You'll probably see a _lot_ of merge conflicts:

```git
Auto-merging nextflow.config
CONFLICT (add/add): Merge conflict in nextflow.config
Auto-merging main.nf
CONFLICT (add/add): Merge conflict in main.nf
Auto-merging environment.yml
CONFLICT (add/add): Merge conflict in environment.yml
Auto-merging docs/usage.md
CONFLICT (add/add): Merge conflict in docs/usage.md
```

Go through each file resolving the merge conflicts carefully.
Many text editors have plugins to help with this task.

It's highly recommended to use a visual tool to help you with this, as it's easy to make mistakes if handling
the merge markers manually when there are so many to deal with.

Once you have resolved all merge conflicts, you can commit the changes and push to the GitHub repo:

```bash
git commit -m "Merged vanilla TEMPLATE branch into main pipeline"
git push origin template_merge
```

The final task is to create a pull request with your changes so that they are included in the upstream repository.
Once your commits are finally merged into the `master` branch, all future automatic template syncing should work.

When new releases of `nf-core/tools` and it's associated template are released, pull-requests will automatically
be created to merge updates in to your pipeline for you.

That's it, you're done! **Congratulations!**

# Fixing a broken `TEMPLATE` branch

If merge conflicts are resolved via the GitHub interface instead of after pulling changes to a fork as described above, the commit history from the `dev` branch will be merged into `TEMPLATE`.
This leads to complex problems in later `TEMPLATE` merges as the later updated `TEMPLATE` branch removes all the pipeline-specific files that were accidentally included in problematic merge, resulting in many (in some cases >100!) of files to resolve conflicts in.

If during one of the automated syncs you see you have an usually high number of changed files you can check whether `dev` was accidentally merged into `TEMPLATE` by looking at how many commits the `TEMPLATE` branch has (should be in the range of 5-15ish, depending on how many template updates your pipeline has had). You can also look at your repository's GitHub Network Graph under the _"Insights"_ tab or even look through the `TEMPLATE` branch commit history to see if there is a commit described as 'Merge branch `dev` into `TEMPLATE`'.

If so, the easiest solution is to start your `TEMPLATE` branch from scratch.

- Clone the main nf-core pipeline repository to your local machine (not your development fork)

  ```bash
  git clone https://github.com/nf-core/<PIPELINE>.git
  cd <pipeline>
  ```

- Next, retrieve the commit hash when the original nf-core template was used to generate pipeline i.e. with `nf-core create`.

  - Assuming you originally started with the nf-core template, you can simply look at your git log from within your repository:

    ```bash
    git checkout TEMPLATE
    git log --reverse
    ```

  - The first commit will then typically represent the original template, with a commit message like `initial template build from nf-core/tools, version 1.9`

- Reset the `TEMPLATE` branch back to this commit, discarding all changes after that

  ```bash
  # Make sure you're definitely have TEMPLATE checked out!
  git reset --hard <hash of first commit after nf-core create>
  ```

- Push this cleaned branch back to the repository - use `--force` to overwrite the history there:

  ```bash
  git push origin TEMPLATE --force
  ```

  - This will then replace the broken `TEMPLATE` branch on GitHub with a nice clean one, which can be viewable by checking the commit history.
  - :warning: You are irreversibly overwriting git history here - make sure that you get the branch names right!

- We can switch back to `dev`, and run `nf-core sync` to do it's magic and get the latest version of the template.

  ```bash
  git checkout dev
  nf-core sync
  ```

- You probably want to now delete your local copy of the pipeline that you checked out from the main nf-core repository.
- On your personal fork of the pipeline you'll want to pull in this fresh template branch:

  ```bash
  cd <path/to/forked/pipeline>
  git branch -D TEMPLATE # Delete the TEMPLATE branch in your fork if you have it
  git remote add upstream git@github.com:nf-core/<PIPELINE>.git  # You might already have this set up?
  git fetch upstream TEMPLATE
  git checkout --track upstream/TEMPLATE
  git push --force
  ```

With this, you're now ready to re-make the pull request from `TEMPLATE` into `dev`, and locally manually resolve conflicts (if required) following the git instructions [above](#merge-template-into-main-branches).
