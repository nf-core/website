---
title: arm64 on Bioconda
subtitle: Getting packages to build for ARM CPUs
parentWeight: 10
---

# Notes on bioconda recipe porting

## About bioconda?

bioconda is a channel of packages (recipes) that are installed via the conda system. Its recipes are conda recipes that are specific to bio-community, and conda-forge is the general purpose base recipes. bioconda and conda have different CI systems and different structures to recipe files, but also enough similarity.

## What we need to do

A successful result is being able to type "conda install {package}" and get the package to install.

## Getting set up

- To use bioconda (required)

  - Install Conda: NB: say yes at the end to get it set things up when you login next time.

    ```
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
    sh ./Miniforge3-Linux-aarch64.sh
    ```

  - Install bioconda
    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

- To develop locally (optional)

  - Add a conda dev environment to do local build / testing of bioconda. Local testing isn't always reliable (often works when the CI doesn't).

    ```
    conda create -n bioconda -c conda-forge -c bioconda bioconda-utils
    conda activate bioconda
    ```

## Enabling packages

First check why your package doesn't install:

```
conda install {name}
```

Either it says not found - so you need to fix that package - or it lists dependencies that are not found and you now have a list to go after.

Most packages we need to fix are in bioconda, but some are in conda-forge. If there is a recipes/{package-name} subdirectory in bioconda-recipes, then this is a bioconda recipe. All bioconda recipes are in the one [github repository](https://github.com/bioconda/bioconda-recipes/)

If the package is not a bioconda recipe, it will be a conda-forge recipe. These are in individual github repositories (feedstocks) of the conda-forge user. They are of the form {package-name}-feedstock.

### Bioconda Recipes

Packages are either generic or not generic. A generic package can have missing non-generic (architecture specific, binary) dependencies and hence not work until the dependencies are made.

If a package is generic, the file recipes/{name}/meta.yaml will have this in the top level build section:

```
build:
..
  noarch: generic
..
```

If a generic package doesn't work - it's a dependency at fault - go fix that. I've recently tried bumping the build number and doing a PR as I think generic packages are locked to versions they were built with - and if those were built before the linux-aarch64 support was added to a dependency, it might not be able to find it - so a rebuild can't harm. I'll update this with the outcome.

If a package is not generic - at some point it compiles native code, or uses binaries that it downloads (and usually for x86 only). Examine the meta.yaml and the adjacent build.sh files.

Normally, to enable it to try to build - add this section to the bottom of the meta.yaml:

```
extra:
 additional-platforms:
   - linux-aarch64
   - osx-arm64
```

you must also bump the build number (eg. add 1 to the existing number), and if the build section is missing a package versioning ('pin') line, you must add one to pass the linter.

```
build:
  number: 2
  noarch: generic
  run_exports:
    - {{ pin_subpackage(name, max_pin='x.x') }}
```

The Pull Request instructions (shown during opening a pull request) will explain when to use "x.x" and when to use "x". In general, I've used 'x' if the major version is > 0.

You are now ready to try it in the CI. Open a Pull Request. Watch the CI and resolve any errors!

- If you need to patch a source file - meta.yaml can handle patches.
- Most recipes have a build.sh - which builds things and can do bash logic to choose different paths for different platforms. CFLAGS and CC etc are set before it runs this script - use those, not the system's gcc. ${PREFIX} is the directory base things are built from - with a bin, lib and include containing all the dependencies specified in meta.yaml
- In the CI output, you'll likely notice that conda packages basically build and run in their own world of lib and bin subdirectories containing everything that they need. You'll see mad looking directories "placehold\_....." - which are just it creating an environment to isolate itself in for build, and another for test.

You can also build and test packages locally to be a bit quicker - it's not 100% reliable to track errors, although the container approach probably is. See [bioconda dev instructions](https://bioconda.github.io/contributor/building-locally.html).

If you want to build outside of the docker container (the containers often have UID and permission errors for me, YMMV)

```
bioconda-utils build --packages {package}
```

### Conda Forge

To fix a conda-forge package - say 'perl-nonsense'

- Head to https://github.com/conda-forge/perl-nonsense-feedstock.
- Check for any open pull requests, one might be an erroring migration to support Arm.
  - If there is one, and if it ran ages ago, the CI logs are probably deleted - so add the comment:
    ```
    @conda-forge-admin please rerender
    ```
    this will trigger a new build so that you can see why it's failing now, if it fails.
- If there isn't an outstanding pull request, you get the system to try migrating for you: edit the [migrations list](https://github.com/conda-forge/conda-forge-pinning-feedstock/blob/main/recipe/migrations/arch_rebuild.txt) to just add the package name - in the right place alphabetically - and get github to make a pull request for you. The bots will open a pull request in the perl-nonsense feedstock and try building the recipe in the next couple of hours or so.
- If successful, the maintainer needs to merge it. After a merge, you still need to wait a couple of days for the binary to be built and made available, and you should be able to do "conda install perl-nonsense".
  - If the maintainer is unresponsive after a few days, try "@conda-forge-admin, please ping team" and if there is no response in a week then I do "@conda-forge/core please help review and merge this PR"
- If the build was unsuccessful, you need to fix it:

```
git clone git@github.com:{gitid}/perl-nonsense-feedstock.git
cd perl-scalar-list-utils-feedstock
git remote add bot  https://github.com/regro-cf-autotick-bot/perl-nonsense-feedstock
git fetch bot
git checkout -b aarch64-fixes bot/bot-pr_arch_[TAB][TAB]
```

- Fix the problem in : edit meta.yaml or conda-forge.yml

```
git commit --all
git push --set-upstream origin ...
```

- Create a new PR from your branch.
- You need to ask the conda-forge bot to do some 'rerendering' which builds a ton of config files / scripts from that recipe - add a comment to your PR of `@conda-forge-admin please rerender`. If you have edited anything, it doesn't harm.
- Note things are very slow at the moment - it tries to build on linux-ppc64le and linux-aarch64 using emulation (!) or Travis's arm fleet which have been erroring. The platforms that are tried is set in the feedstock's conda-forge.yml file - if ppc64le is failing, you can remove the entry and just fix the linux-aarch64 one.
