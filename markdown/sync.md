
# Synchronisation
To keep all the nf-core pipelines up-to-date with the latest version of the community standards, we have implemented a synchronisation tool.
This ensures that updates to the community standards are propagated to all the nf-core pipelines.

If a pipeline is created using the `nf-core create` command, everything is automagically set up for the synchronisation to work.
However, if a pipeline was not created using this command, this documentation will tell you how to set it up properly.

## How it works
In brief, when the synchronisation tool is triggered it compares every pipeline against the new standards and opens pull requests when mismatches are found.
These pull requests needs to be manually dealt with by the pipeline maintainers.

Behind the scenes, the synchronisation tool fetches the variables needed for a pipeline and uses this to trigger a <nobr>`nf-core create --no-git`</nobr> command.
The result from this command is then compared against what is stored in a special branch of the pipeline called the `TEMPLATE` branch.
If any changes are detected, a pull request is created.

For this to work in practice, the `TEMPLATE` branch needs to have a shared git history with the `master` branch of the pipeline.
The <nobr>`nf-core create`</nobr> command arranges this by enforcing a first commit to the `master` branch before any development has taken place.
The next section deals with the case when the pipeline *was not* created by the `nf-core create` command.

## <a name="setup"></a> Set up your pipeline for automatic syncing
This section describes how to set up a correct TEMPLATE branch in the case your pipeline was not created with a TEMPLATE branch from the beginning. If you created a pipeline with the `nf-core create` command, you should be all ready to go and can skip this step. Otherwise proceed with caution. It is probably a good idea to make sure you have all your local changes pushed to github and you could even make a local backup clone of your repository before proceeding.

You should also consider the option to restart your pipeline project by running the `nf-core create` command and simply copy in the modifications you need into the newly created pipeline.

### Step-by-step procedure
Create TEMPLATE branch:
```bash
cd pipeline_root_dir
git checkout --orphan TEMPLATE
```

Now remove all the files of your pipeline to be able to have a completely empty branch.
```bash
git rm -rf '*'
```
Make sure your branch is completely empty by checking the status of `git status`, which should be completely empty:
```bash
$ git status
On branch TEMPLATE

No commits yet

nothing to commit (create/copy files and use "git add" to track)
```

Regenerate your pipeline from scratch using the most recent template:
*Make sure you are within your pipeline root directory before running these commands.*
```bash
nf-core create --no-git -n 'YOURPIPELINENAME' -d 'YOUR ACTUAL PIPELINE DESCRIPTION'
```

This creates a new directory `YOURPIPELINENAME` with the template pipeline files in it.
Now move these files into your root git directory:
```bash
mv nf-core-YOURPIPELINENAME/* .
mv nf-core-YOURPIPELINENAME/.[!.]* .
rm -rf nf-core-YOURPIPELINENAME
```
Now make sure the newly created files are in the correct place. It should look similar to this:
```
$ git status
On branch TEMPLATE

No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	.gitattributes
	.github/
	.gitignore
	.travis.yml
	CHANGELOG.md
	CODE_OF_CONDUCT.md
	Dockerfile
	LICENSE
	README.md
	Singularity
	assets/
	bin/
	conf/
	docs/
	environment.yml
	main.nf
	nextflow.config

nothing added to commit but untracked files present (use "git add" to track)
```
If it all looks good, then commit these files:
```
git add .
git commit -m “Initial template commit”
```
For the nf-core bot to be able to access your TEMPLATE branch, you need to push it to the upstream repository.
```
git push upstream TEMPLATE
```

### Merge TEMPLATE into master
The only remaining step is unfortunately a rather tedious one.
You have to merge the TEMPLATE branch into your master (or development) branch and manually resolve all merge conflicts.
```
git checkout master
git merge TEMPLATE --allow-unrelated-histories
git status
```
Make sure to go through all the conflicts and resolve them so that your pipeline is properly up to date with the latest template.
Finally, commit and upload your changes to your github repo.
```
git add .
git commit -m "Merged vanilla TEMPLATE branch into master"
git push origin master
```
The final task is to create a pull request with your changes so that they are included in the upstream repository.

## *Well done!*
