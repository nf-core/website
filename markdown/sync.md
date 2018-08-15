
# Advanced users - Adherence to the nf-core template
When a new pipeline is created using the command `nf-core create`, quite a few things happen behind the scene.
One of those things is to setup a git branch called TEMPLATE for your pipeline.
This is branch is used to track the changes made to the central nf-core template used by `nf-core create`.
The reason to keep a git branch for this on each pipeline is to enable updates to the central template to be merged into individual pipelines.

## Set up automatic adherence to the nf-core template
This section describes how to set up a correct TEMPLATE branch in the case your pipeline was not created with a TEMPLATE branch from the beginning. If you created a pipeline with the `nf-core create` command, you should be all set up and should skip this step. Otherwise proceed with caution. It is probably a good idea to make sure you have all your local changes pushed to github and you could even make a local backup clone of your repository before proceeding.

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
Finally, commit and upload your changes to your github repo and create a pull request.
```
git add .
git commit -m "Merged vanilla TEMPLATE branch into master"
git push origin master
```

*Well done!*
