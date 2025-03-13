# nf_core.lint

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to
the nf-core community guidelines.

### _`class{:python}`_`nf_core.lint.PipelineLint(path){:python}`

Object to hold linting information and results.
All objects attributes are set, after the [`PipelineLint.lint_pipeline()`](#nf_core.lint.PipelineLint.lint_pipeline) function was called.

- **Parameters:**
  **path** (_str_) – The path to the nf-core pipeline directory.

#### `conda_config{:python}`

The parsed conda configuration file content (environment.yml).

- **Type:**
  dict

#### `conda_package_info{:python}`

The conda package(s) information, based on the API requests to Anaconda cloud.

- **Type:**
  dict

#### `config{:python}`

The Nextflow pipeline configuration file content.

- **Type:**
  dict

#### `dockerfile{:python}`

A list of lines (str) from the parsed Dockerfile.

- **Type:**
  list

#### `failed{:python}`

A list of tuples of the form: (<error no>, <reason>)

- **Type:**
  list

#### `files{:python}`

A list of files found during the linting process.

- **Type:**
  list

#### `minNextflowVersion{:python}`

The minimum required Nextflow version to run the pipeline.

- **Type:**
  str

#### `passed{:python}`

A list of tuples of the form: (<passed no>, <reason>)

- **Type:**
  list

#### `path{:python}`

Path to the pipeline directory.

- **Type:**
  str

#### `pipeline_name{:python}`

The pipeline name, without the nf-core tag, for example hlatyping.

- **Type:**
  str

#### `release_mode{:python}`

True, if you the to linting was run in release mode, False else.

- **Type:**
  bool

#### `warned{:python}`

A list of tuples of the form: (<warned no>, <reason>)

- **Type:**
  list

**Attribute specifications**

Some of the more complex attributes of a PipelineLint object.

- conda_config:
  ```default
  # Example
   {
      'name': 'nf-core-hlatyping',
      'channels': ['bioconda', 'conda-forge'],
      'dependencies': ['optitype=1.3.2', 'yara=0.9.6']
   }
  ```
- conda_package_info:
  ```default
  # See https://api.anaconda.org/package/bioconda/bioconda-utils as an example.
   {
      <package>: <API JSON repsonse object>
   }
  ```
- config: Produced by calling Nextflow with `nextflow config -flat <workflow dir>`. Here is an example from
  : the [nf-core/hlatyping](https://github.com/nf-core/hlatyping) pipeline:
  ```default
  process.container = 'nfcore/hlatyping:1.1.1'
  params.help = false
  params.outdir = './results'
  params.bam = false
  params.singleEnd = false
  params.seqtype = 'dna'
  params.solver = 'glpk'
  params.igenomes_base = './iGenomes'
  params.clusterOptions = false
  ...
  ```

#### `check_anaconda_package(dep){:python}`

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

- **Parameters:**
  **dep** (_str_) – A conda package name.
- **Raises:**
  **A ValueError**\*\*,\*\* **if the package name can not be resolved.** –

#### `check_ci_config(){:python}`

Checks that the Travis or Circle CI YAML config is valid.

Makes sure that `nf-core lint` runs in travis tests and that
tests run with the required nextflow version.

#### `check_conda_dockerfile(){:python}`

Checks the Docker build file.

Checks that:
: \* a name is given and is consistent with the pipeline name

- dependency versions are pinned
- dependency versions are the latest available

#### `check_conda_env_yaml(){:python}`

Checks that the conda environment file is valid.

Checks that:
: \* a name is given and is consistent with the pipeline name

- check that dependency versions are pinned
- dependency versions are the latest available

#### `check_docker(){:python}`

Checks that Dockerfile contains the string `FROM`.

#### `check_files_exist(){:python}`

Checks a given pipeline directory for required files.

Iterates through the pipeline’s directory content and checkmarks files
for presence.
Files that **must** be present:

```default
'nextflow.config',
'Dockerfile',
['.travis.yml', '.circle.yml'],
['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling
'README.md',
'CHANGELOG.md',
'docs/README.md',
'docs/output.md',
'docs/usage.md'
```

Files that _should_ be present:

```default
'main.nf',
'environment.yml',
'conf/base.config'
```

- **Raises:**
  **An AssertionError if neither nextflow.config** **or** **main.nf found.** –

#### `check_licence(){:python}`

Checks licence file is MIT.

Currently the checkpoints are:
: \* licence file must be long enough (4 or more lines)

- licence contains the string _without restriction_
- licence doesn’t have any placeholder variables

#### `check_nextflow_config(){:python}`

Checks a given pipeline for required config variables.

Uses `nextflow config -flat` to parse pipeline `nextflow.config`
and print all config variables.
NB: Does NOT parse contents of main.nf / nextflow script

#### `check_pip_package(dep){:python}`

Query PyPi package information.

Sends a HTTP GET request to the PyPi remote API.

- **Parameters:**
  **dep** (_str_) – A PyPi package name.
- **Raises:**
  **A ValueError**\*\*,\*\* **if the package name can not be resolved** **or** **the connection timed out.** –

#### `check_pipeline_todos(){:python}`

Go through all template files looking for the string ‘TODO nf-core:’

#### `check_readme(){:python}`

Checks the repository README file for errors.

Currently just checks the badges at the top of the README.

#### `check_version_consistency(){:python}`

Checks container tags versions.

Runs on `process.container`, `process.container` and `$TRAVIS_TAG` (each only if set).

Checks that:
: \* the container has a tag

- the version numbers are numeric
- the version numbers are the same as one-another

#### `lint_pipeline(release_mode=False){:python}`

Main linting function.

Takes the pipeline directory as the primary input and iterates through
the different linting checks in order. Collects any warnings or errors
and returns summary at completion. Raises an exception if there is a
critical error that makes the rest of the tests pointless (eg. no
pipeline script). Results from this function are printed by the main script.

- **Parameters:**
  **release_mode** (_boolean_) – Activates the release mode, which checks for
  consistent version tags of containers. Default is False.
- **Returns:**
  Summary of test result messages structured as follows:
  ```default
  {
      'pass': [
          ( test-id (int), message (string) ),
          ( test-id (int), message (string) )
      ],
      'warn': [(id, msg)],
      'fail': [(id, msg)],
  }
  ```
- **Return type:**
  dict
- **Raises:**
  **If a critical problem is found**\*\*,\*\* **an AssertionError is raised.** –

### `nf_core.lint.run_linting(pipeline_dir, release_mode=False){:python}`

Runs all nf-core linting checks on a given Nextflow pipeline project
in either release mode or normal mode (default). Returns an object
of type [`PipelineLint`](#nf_core.lint.PipelineLint) after finished.

- **Parameters:**
  - **pipeline_dir** (_str_) – The path to the Nextflow pipeline root directory
  - **release_mode** (_bool_) – Set this to True, if the linting should be run in the release mode.
    See [`PipelineLint`](#nf_core.lint.PipelineLint) for more information.
- **Returns:**
  An object of type [`PipelineLint`](#nf_core.lint.PipelineLint) that contains all the linting results.
