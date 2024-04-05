---
title: Fixing a broken TEMPLATE branch
description: How to fix a broken TEMPLATE branch in your nf-core pipeline repository
---

If merge conflicts are resolved via the GitHub interface instead of after pulling changes to a fork as described [in the `sync` documentation](/docs/contributing/sync), the commit history from the `dev` branch will be merged into `TEMPLATE`.
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

  :::warning
  You are irreversibly overwriting git history here - make sure that you get the branch names right!
  :::

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

With this, you're now ready to re-make the pull request from `TEMPLATE` into `dev`, and locally manually resolve conflicts (if required) following the git instructions [in the `sync` documentation](/docs/contributing/sync#merge-template-into-main-branches).
