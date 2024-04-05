---
title: Using GitHub to contribute to nf-core
subtitle: A guide to contributing to nf-core using GitHub
---

# Fork the nf-core/modules repository and branch

The first step, to contribute a module to the community repository is to fork \*nf-core modules into your own account or organisation. To do this, you should click on the top-right of the nf-core modules repository, and choose "fork" as shown in the figure below.

![fork](/images/contributing/dsl2_modules_tutorial/dsl2-mod_01_fork.png)

You then choose the account or organisation you want to fork the repository into. Once forked, you can commit all changes you need into the new repository.

In order to create a new module, it is best to branch the code into a recognisable branch. You can do this in two ways.

- You can create a new branch locally, on the terminal, using the following command:

  - ```bash
    git checkout -b newmodule

    ## alternatively: git switch -c newmodule
    ```

  - The branch will be synchronised with your remote once you push the first new commit.

- You can use the GitHub interface

  - To do this, you can select the dropdown menu on the top-left of your repository code, write the name of the new branch and choose to create it as shown below:

    ![branch](/images/contributing/dsl2_modules_tutorial/dsl2-mod_02_new_branch.png)

  - You will then sync this locally (ideally, you clone the forked repository on your working environment to edit code more comfortably)
