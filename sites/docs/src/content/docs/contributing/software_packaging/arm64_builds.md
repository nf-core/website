---
title: arm64 on Bioconda
subtitle: Getting packages to build for ARM CPUs
shortTitle: arm64 builds
parentWeight: 10
---

## Bioconda and Conda-forge

Bioconda is a channel of software packages (recipes) that are installed via the conda system.
Its recipes are conda recipes that are specific to bio-community, and conda-forge is the general purpose base recipes.
Bioconda and conda-forge have different CI systems and different structures to recipe files, but they also have a lot in common.

This page has documentation about taking an existing package that's already on bioconda / conda-forge and making it work not just on `linux/amd64` (intel chips) but also `linux/arm64` chips (ARM chips, like AWS Graviton).

## What we need to do

A successful result is being able to type `conda install {package}` and get the package to install.

## Getting set up

To use Bioconda _(required)_

1. Install Conda:

   ```bash
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
   sh ./Miniforge3-Linux-aarch64.sh
   ```

   Input **Yes** at the end of the install for conda to load when you open a new terminal.

2. Install Bioconda:

   ```bash
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```

To develop locally _(if running on an arm64 machine - optional)_

1. Add a conda dev environment to do local build / testing of Bioconda:

   ```bash
   conda create -n bioconda-builds -c conda-forge -c bioconda bioconda-utils
   conda activate bioconda-builds
   ```

:::note
Local testing isn't always reliable. It often works when the CI doesn't.
:::

## Enabling packages

Packages can be tested locally or by using Wave.

### Testing locally

If you are running on an ARM machine, you can first check why your package doesn't install:

```bash
conda install {name}
```

### Testing using Wave

If not, you can try to build using [Wave](https://seqera.io/wave/) - either by requesting a container via [Seqera Containers](https://seqera.io/containers/) (remember to select `linux/arm64` in the settings) or by using the [Wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda {name} --platform linux/arm64 --freeze --await
```

The `--freeze` flag tells the CLI to store the generated images for everyone to use on Seqera Containers.
The `--await` flag tells the CLI to keep running until the build is complete.
If the build fails when using `--await`, you'll get an error that looks something like this:

```
Container provisioning did not complete successfully
- Reason: Container build did not complete successfully
- Find out more here: https://wave.seqera.io/view/builds/bd-xxxxxxx
```

If you see this error, follow the link to go to the _build details_ web page. It includes the full conda output, as if you'd run it locally.

### Interpreting output

If the package doesn't build, there are usually two types of error:

1. It says "not found" - you need to fix that package
2. It lists dependencies that are not found - you now have a list of dependencies to investigate.

Most packages that need to be fixed are in Bioconda. However, some are in conda-forge.
If there is a `recipes/{package-name}` subdirectory in the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/) GitHub repo, then this is a Bioconda recipe.

If the package is not a Bioconda recipe, it will be a conda-forge recipe.
These are in individual GitHub repositories (feedstocks) under the [conda-forge GitHub organisation](https://github.com/conda-forge/).
Each repository is named in the form `{package-name}-feedstock`.

## Bioconda Recipes

Bioconda recipes can be generic or non-generic.

### Generic recipes

A generic package can have missing non-generic (architecture-specific, binary) dependencies and, hence, not work until the dependencies are made.

If a package is generic, the file `recipes/{name}/meta.yaml` will have this in the top level build section:

```yaml
build:
..
  noarch: generic
..
```

If a generic package doesn't work - a dependency is at fault and will require fixing.

This can be as simple as bumping the build number and doing a PR, as generic packages are locked to versions they were built with.
If those were built before the bioconda `linux-aarch64` support was added to a dependency, it might not be findable - so a rebuild won't harm it.

### Not-generic

If a package is not generic, at some point, it compiles native code or uses binaries that it downloads (and usually for `x86` only).
Examine the `meta.yaml` and the adjacent `build.sh` files.

Normally, to enable it to try to build, add this section to the bottom of the `meta.yaml`:

```yaml
extra:
  additional-platforms:
    - linux-aarch64
    - osx-arm64
```

**You must also bump the build number** (For example, add `1` to the existing number).
If the build section is missing a package versioning (`pin`) line, you must add one to pass the linter.

```yaml
build:
  number: 2
  noarch: generic
  run_exports:
    - { { pin_subpackage(name, max_pin='x.x') } }
```

The Pull Request instructions (shown during opening a pull request) will explain when to use `'x.x'` and when to use `'x'`.
In general, I've used `'x'` if the major version is > 0.

### Testing in the CI

You are now ready to try it in the CI. Open a Pull Request. Watch the CI and resolve any errors!

- If you need to patch a source file - `meta.yaml` can handle patches.
- Most recipes have a `build.sh` which builds things and can do bash logic to choose different paths for different platforms.
  `CFLAGS` and `CC` etc are set before it runs this script - use those, not the system's `gcc`.
  `${PREFIX}` is the directory base things are built from - with a `bin`, `lib` and `include` containing all the dependencies specified in `meta.yaml`
- In the CI output, you'll likely notice that conda packages basically build and run in their own world of `lib` and `bin` subdirectories containing everything that they need.
  You'll see mad looking directories `"placehold_....."` - which are just it creating an environment to isolate itself in for build, and another for test.

### Testing builds locally

You can also build and test packages locally to be a bit quicker.
It's not 100% reliable to track errors, although the container approach probably is.
See [Bioconda dev instructions](https://bioconda.github.io/contributor/building-locally.html).

If you want to build outside of the docker container (the containers often have UID and permission errors for me, YMMV), the command is:

```bash
bioconda-utils build --packages {package}
```

## Conda Forge

### Quick fixes

To fix a conda-forge package, for example one called `perl-nonsense`:

1. Go to https://github.com/conda-forge/perl-nonsense-feedstock
2. Check for any open pull requests, one might be an erroring migration to support Arm.

- If there is one, and if it ran ages ago, the CI logs are probably deleted - so add the comment:
  ```
  @conda-forge-admin please rerender
  ```
  This will trigger a new build so that you can see why it's failing now, if it fails.

3. If there isn't an outstanding pull request, you get the system to try migrating for you:

- Edit the [migrations list](https://github.com/conda-forge/conda-forge-pinning-feedstock/blob/main/recipe/migrations/arch_rebuild.txt) to just add the package name - in the right place alphabetically - and get github to make a pull request for you.
- The bots will open a pull request in the `perl-nonsense` feedstock and try building the recipe in the next couple of hours or so.

4. If successful, the maintainer needs to merge it.

- After a merge, you still need to wait a couple of days for the binary to be built and made available, and you should be able to do:
  ```
  conda install perl-nonsense
  ```
- If the maintainer is unresponsive after a few days, try adding a comment with the phrase:
  ```
  @conda-forge-admin, please ping team
  ```
- If there is still no response, try adding a comment with the phrase:
  ```
  @conda-forge/core please help review and merge this PR
  ```

5. If the build was unsuccessful, you need to fix it

### Manual fixes

To fix a recipe, you need to fork the repo and clone it locally:

```bash
git clone git@github.com:{your-github-username}/perl-nonsense-feedstock.git
cd perl-nonsense-feedstock
git remote add bot https://github.com/regro-cf-autotick-bot/perl-nonsense-feedstock
git fetch bot
git checkout -b aarch64-fixes bot/bot-pr_arch_[TAB][TAB]
```

Then fix the problem: edit `meta.yaml` or `conda-forge.yml` and commit / push:

```bash
git commit --all
git push --set-upstream origin aarch64-fixes
```

Next, create a new PR from your branch.

You need to ask the conda-forge bot to do some 'rerendering' which builds a ton of config files / scripts from that recipe:

- Add a comment to your PR of `@conda-forge-admin please rerender`. If you have edited anything, it doesn't harm.
- Note things are very slow at the moment
  - It tries to build on `linux-ppc64le` and `linux-aarch64` using emulation (!) or Travis's arm fleet which have been erroring.
  - The platforms that are tried is set in the feedstock's `conda-forge.yml` file - if `ppc64le` is failing, you can remove the entry and just fix the `linux-aarch64` one.

Done!
