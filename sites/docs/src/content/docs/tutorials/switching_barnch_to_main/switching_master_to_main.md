---
title: "Transitioning from `master` to `main` branches in pipelines"
subtitle: Guide on how to switch a pipeline from the default branch `master` to `main`
weight: 20
---

This tutorial will guide you through changing the pipeline default branch form `master` to `main`.

## Prequisites

- Admin rights on the GitHub repository where the change is taking place
- nf-core tools
  - It is recommended that the version match the version of the template in your pipeline

> [!WARNING]
> Make sure you are on using same version of nf-core/tools as the template version of your branch.
> This will avoid adding other changes not needed for the branch name switch.

## Instructions

1. Open your repository on GitHub (make sure you're on `master`)
2. Go to the `branches` page (to the right of the `branches` drop down where it says 'master')

![image](https://hackmd.io/_uploads/rkc8GXJWgx.png)

3. Under the default branch section of the page, press the triple dot menu on the `master` branch row, and press 'Rename branch'

:::note
You will need permissions to change this setting. If you don't have permission, ask the @core-team to do it for you.
:::

![image](https://hackmd.io/_uploads/Bkpdz7kWlg.png)

4.  Rename `master` to `main`
5.  Press the 'learn more' text, and copy and paste the displayed instructions somewhere safe.

> [NOTE]
> These instructions will be usefur for all pipeline contributors.

![image](https://hackmd.io/_uploads/SkiIQXyWlx.png)

6. Press 'Rename branch'
7. Go back to 'Code' tab to verify you are now on `main`
8. In your local IDE (e.g. VSCode) make sure you're on the `dev` branch
9. Run `git pull` to ensure you have the `main` branch locally
10. Check your `git` config to see what is your current default branch (it should report `master`):

    ```bash
    git config --global init.defaultBranch
    ```

11. Change it to `main` with :

    ```bash
    git config --global init.defaultBranch main
    ```

12. Check your `git` config again to check it change (it should now report `main`):

    ```bash
    git config --global init.defaultBranch
    ```

13. Still on `dev`, run `nf-core pipelines sync`
14. Change to a new branch `git switch -c default-branch-change`
15. Follow the merge TEMPLATE instructions, i.e. `git merge TEMPLATE`

> [NOTE]
> If you don't want any other template changes, make sure to use `nf-core/tools` version that matches the template version in your pipeline

18. Resolve merge conflicts

> [TIP]
> If it's the ROcrate file, you can accept all incoming chagnes

19. Check that all references of your pipelines' `msater` is now `main` using a global repository search in your IDE

> [!WARNING]
> Make sure not to modify references of master in links to other repositories!
> If in doubt, ask on the nf-core slack!

20. Run `nf-core pipelines lint` to check you didn't break anything
21. Commit and merge `git add -am 'Change default branch'
22. Push the changes `git push`
23. On GitHub make a new PR against `dev`
24. Review PR to check all links and `nextflow.config` manifest say `main`
25. Review and merge :tada:

## Post-change instructions

1. Inform all your collaborators (on slack, etc.) that the default branch has changed, and that they should update their clones and forks, i.e.:

   ```bash
   git branch -m master main
   git fetch origin
   git branch -u origin/main main
   git remote set-head origin -a
   ```

2. Update your own clones and forks
