---
title: Pipeline Development Environment
subtitle: Where and how we develop our pipelines
---

## GitHub usage

### Repository organisation

The default branch of our repository is called `main`, following the new GitHub convention, as opposite to `master` in [nf-core/modules](https://github.com/nf-core/modules).

The `main` and `dev` branches are "protected" and can only be modified with pull-requests.

### Access

- All our pipeline repositories are public: everyone has read access.
- All sanger-tol Nextflow developers have write access to all Nextflow repositories.
  - Direct push to `main` and `dev` is not allowed, since the two branches are protected.
- Within each team (Tree of Life Assembly and Informatics Infrastructure), everyone has ["maintain" access](https://docs.github.com/en/organizations/managing-user-access-to-your-organizations-repositories/managing-repository-roles/repository-roles-for-an-organization#permissions-for-each-role) to the team repositories.
- A selection of people from the Informatics Infrastructure team have ["admin" access](https://docs.github.com/en/organizations/managing-user-access-to-your-organizations-repositories/managing-repository-roles/repository-roles-for-an-organization#permissions-for-each-role) to all Nextflow repositories.
- @muffato and @mcshane are ["owners"](https://docs.github.com/en/organizations/managing-peoples-access-to-your-organization-with-roles/roles-in-an-organization#permissions-for-organization-roles).
- Additional contributors can be added, if sponsored by a repository admin.

## Farm environment

To develop and run pipelines on our LSF compute farms, first make sure you have set up the "module" and "Conda" environments as per https://ssg-confluence.internal.sanger.ac.uk/display/TOL/Farm+environment

Modules are used to expose Nextflow alone

```
$ module load nextflow/23.10.0-5889
$ nextflow -version

      N E X T F L O W
      version 23.10.0 build 5889
      created 15-10-2023 15:07 UTC (16:07 BST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

Conda environments are used to expose the nf-core command-line tool, together with a Nextflow and some development helpers like `prettier` or `nf-test`.

```
$ conda activate nf-core_2.11
$ nextflow -version

      N E X T F L O W
      version 23.04.1 build 5866
      created 15-04-2023 06:51 UTC (07:51 BST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io

$ nf-core --version

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.11.1 - https://nf-co.re
    There is a new version of nf-core/tools available! (2.13.1)


nf-core, version 2.11.1
$ nf-test version
🚀 nf-test 0.8.2
https://code.askimed.com/nf-test
(c) 2021 - 2023 Lukas Forer and Sebastian Schoenherr

Nextflow Runtime:
    >
      N E X T F L O W
      version 23.04.1 build 5866
      created 15-04-2023 06:51 UTC (07:51 BST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io

$ prettier --version
3.1.1
```

Java programs are not really allowed to run on the head nodes so you will have to submit all your Nextflow commands on LSF.
When developing a pipeline, you may want to do all your development from an interactive job.

In our experience, `nextflow run` commands need 1 CPU and 6 GB RAM for themselves, but of course more if you want to use the "local" executor.
For all sanger-tol and nf-core pipelines, we recommend using the ["sanger" profile](https://github.com/nf-core/configs/blob/master/conf/sanger.config)
which automatically activates LSF job submission using the appropriate queues based on each job's parameters.
When asking Nextflow to submit jobs to LSF, please submit Nextflow itself to the `oversubscribed` queue so that it doesn't take unnecessary compute resources.
