# nf-core tools

To help you get started with nf-core, we've put together a companion package, written in Python.

## Installation

> NB: This is not correct - nf-core hasn't been released yet. Soon!

```
pip install nf-core
```

Alternatively, you can install the development version directly from GitHub:

```
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git
```

## <a name="create"></a>Create a pipeline from scratch 
To create a pipeline from scrath, use the command `nf-core create`.

You should see an output like this:

```
$ nf-core create

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'
    
Usage: nf-core create [OPTIONS]

  Create a new pipeline using the nf-core template

Options:
  -n, --name TEXT         The name of your new pipeline  [required]
  -d, --description TEXT  A short description of your pipeline  [required]
  --new-version TEXT      The initial version number to use
  --no-git                Do not initialise pipeline as new git repository
  -f, --force             Overwrite output directory if it already exists
  -o, --outdir TEXT       Output directory for new pipeline (default: pipeline
                          name)
  --help                  Show this message and exit.
```

If you have already initialized a Git repository, you can pass the option flag `--no-git`. Otherwise,
nf-core tools will create a Git repository via `git init`. You still have to add remote tracking by yourself though.

One very important setup is the creation of a branch `TEMPLATE`, right after the initialisation
from the nf-core template. We will use this branch to keep pipelines in sync with future changes
in the template. 

Please follow the [template branch](template_branch.md) description.

## Local development and testing

If you have created your pipeline with the `nf-core create` command, there is already a Dockerfile created automatically. The container receipe will install `conda` as package manager, and you can define software
dependencies by adding packages to the `environment.yml` configuration file in the root directory of your nf-core pipeline.

If you now want to **test** the Docker container **locally**, you can do so with the `Docker` command: 

```
$ docker build -t nf-core/<pipeline-name> .
```

Be sure that you have [Docker](https://www.docker.com/) installed locally, and that the Docker daemon is running.

Once the container is build, it is ready to use in your local Docker container registry. You can now run your `nf-core` pipeline with

```
$ nextflow run path/to/my/main.nf -profile test -with-docker
```

**Note:** Make sure, that you have tagged your container locally with the same name as it is in the `nextflow.config` assignment of `params.container`.


## Listing pipelines
To see available pipelines, use the command `nf-core list`

You should get output that looks something like this:

```
$ nf-core list

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'


Name               Version    Published    Last Pulled    Default local is latest release?
-----------------  ---------  -----------  -------------  ----------------------------------
nf-core/methylseq  1.0        2 days ago   2 days ago     No
nf-core/ExoSeq     dev        -            -              No
nf-core/RNAseq     dev        -            -              No
```

## Lint tests
If you're writing a new pipeline or modifying an existing one, you should run the nf-core lint tests locally before pushing any changes. These tests are run automatically using the Travis CI testing and must pass before pull-requests can be merged.

To run the tests, us the command `nf-core lint [path]`, where `[path]` is the path to your pipeline code. For example, if your current working directory is the pipeline:

```
$ nf-core lint .

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

Running pipeline tests  [####################################]  100%  None

INFO: ===========
 LINTING RESULTS
=================
  69 tests passed   0 tests had warnings   0 tests failed
```

## Bumping version numbers
Changing the version number for a pipeline is a little tedious and mistakes are easily made. To make this a little faster and simpler, there is a command `nf-core release [path] [version]`. For example:

This command runs the lint tests first, to ensure that the existing code meets the expected criteria, before trying to replace exiting variables.

```
$ nf-core release . 1.4

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'


INFO: Running nf-core lint tests
Running pipeline tests  [####################################]  100%  None

INFO: ===========
 LINTING RESULTS
=================
  52 tests passed   0 tests had warnings   0 tests failed

INFO: Changing version number:
  Current version number is '0.1.0'
  New version number will be '1.4'

INFO: Updating version in nextflow.config
 - version = "0.1.0"
 + version = '1.4'

INFO: Updating version in nextflow.config
 - container = 'nfcore/example:0.1.0'
 + container = 'nfcore/example:1.4'

INFO: Updating version in environment.yml
 - name: nfcore-example-0.1.0
 + name: nfcore-example-1.4

INFO: Updating version in Dockerfile
 - ENV PATH /opt/conda/envs/nfcore-example-0.1.0/bin:$PATH
 + ENV PATH /opt/conda/envs/nfcore-example-1.4/bin:$PATH
```



