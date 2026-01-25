---
title: Merging automated PRs
subtitle: How to merge automated PRs after a new nf-core/tools release
shortTitle: Merging automated PRs
---

When nf-core releases a new tools version, each pipeline receives an automated pull request to merge template changes into the pipeline.

If the PR has no merge conflicts, review the changes and merge it into the `dev` branch directly.

If the PR contains merge conflicts, you need to resolve them manually.
You can either work on the branch created for the template sync through the GitHub web interface, or pull the updates to `TEMPLATE` on your local fork.

For non-trivial merges, work on your fork locally.
Comment on the automated PR to notify others that you are resolving the conflicts.
See [Resolve major conflicts](#resolve-major-conflicts) for guidance.

## Resolve minor conflicts

Use this approach for simple conflicts resolvable through the GitHub web interface.

1. Navigate to the **Pull Requests** tab of your repository.

1. Open the PR typically named `Important! Template update for nf-core/tools v<version>`. This PR comes from a branch named `nf-core-template-merge-<version>`, which is a modifiable copy of `TEMPLATE`.

1. Resolve the conflicts at the bottom of the PR page.

1. Request reviews from the nf-core community once tests pass.

## Resolve major conflicts

Use this approach for large conflicts that cannot be resolved through the GitHub web interface.
Local resolution allows you to test changes on your machine before committing.

### Pull the changes to your fork

1. Navigate to the directory where you have checked out your fork of the pipeline repository.

1. Add the nf-core repository as a git remote called `upstream`:

   ```bash
   git remote add upstream https://github.com/nf-core/<pipeline_name>.git
   ```

1. Check out a new branch for these changes:

   ```bash
   git checkout -b merging-template-updates
   ```

1. Pull the `TEMPLATE` branch from the `upstream` repository:

   ```bash
   git pull upstream TEMPLATE
   ```

### Resolve merge conflicts

When you pull the `TEMPLATE` branch, you may see log messages indicating merge conflicts:

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

Check the current status to see files with merge conflicts under **Unmerged paths**:

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

Go through each file to resolve every merge conflict.
Most code editors provide tools to help with this.
For example, [VS Code](https://code.visualstudio.com/docs/editor/versioncontrol#_merge-conflicts) has built-in merge conflict support.

:::note
Use the version from the `TEMPLATE` branch in most cases.
However, some template code may need customisation for your pipeline.
You may need to manually combine both versions into one code block.
:::

If you need assistance, ask for help on the nf-core Slack.

### Push the resolved changes to your fork

Once you resolve all merge conflicts and stage all files, commit and push the changes:

```bash
git commit -m "Merged changes from nf-core template"
git push --set-upstream origin merging-template-updates
```

### Create a pull request

1. Create a pull request from your fork to the main nf-core repository for the pipeline.

1. Request reviews and merge as usual.

The commit history on your PR will include a commit by the `@nf-core-bot` user with the same commit hash as the automated `TEMPLATE` PR.
Once you merge your fork, the automated PR will show as merged and close automatically.
