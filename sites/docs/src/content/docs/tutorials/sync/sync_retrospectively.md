---
title: Setting up a pipeline for syncing retrospectively
subtitle: How to set up a correct TEMPLATE branch in your nf-core pipeline repository
---

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
