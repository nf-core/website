# nf_core.lint

Linting policy for nf-core pipeline projects.

Tests Nextflow-based pipelines to check that they adhere to
the nf-core community guidelines.

### _class_ nf_core.lint.PipelineLint(path)

Object to hold linting information and results.
All objects attributes are set, after the [`PipelineLint.lint_pipeline()`](#nf_core.lint.PipelineLint.lint_pipeline) function was called.

- **Parameters:**
  **path** (_str_) – The path to the nf-core pipeline directory.

#### conda_config

The parsed conda configuration file content (environment.yml).

- **Type:**
  dict

#### conda_package_info

The conda package(s) information, based on the API requests to Anaconda cloud.

- **Type:**
  dict

#### config

The Nextflow pipeline configuration file content.

- **Type:**
  dict

#### dockerfile

A list of lines (str) from the parsed Dockerfile.

- **Type:**
  list

#### failed

A list of tuples of the form: (<error no>, <reason>)

- **Type:**
  list

#### files

A list of files found during the linting process.

- **Type:**
  list

#### minNextflowVersion

The minimum required Nextflow version to run the pipeline.

- **Type:**
  str

#### passed

A list of tuples of the form: (<passed no>, <reason>)

- **Type:**
  list

#### path

Path to the pipeline directory.

- **Type:**
  str

#### pipeline_name

The pipeline name, without the nf-core tag, for example hlatyping.

- **Type:**
  str

#### release_mode

True, if you the to linting was run in release mode, False else.

- **Type:**
  bool

#### warned

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
  params.single_end = false
  params.seqtype = 'dna'
  params.solver = 'glpk'
  params.igenomes_base = './iGenomes'
  params.clusterOptions = false
  ...
  ```

#### check_actions_awsfulltest()

Checks the GitHub Actions awsfulltest is valid.

Makes sure it is triggered only on `release`.

#### check_actions_awstest()

Checks the GitHub Actions awstest is valid.

Makes sure it is triggered only on `push` to `master`.

#### check_actions_branch_protection()

Checks that the GitHub Actions branch protection workflow is valid.

Makes sure PRs can only come from nf-core dev or ‘patch’ of a fork.

#### check_actions_ci()

Checks that the GitHub Actions CI workflow is valid

Makes sure tests run with the required nextflow version.

#### check_actions_lint()

Checks that the GitHub Actions lint workflow is valid

Makes sure `nf-core lint` and `markdownlint` runs.

#### check_anaconda_package(dep)

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

- **Parameters:**
  **dep** (_str_) – A conda package name.
- **Raises:**
  **A ValueError\*\***,\*\* **if the package name can not be resolved.** –

#### check_conda_dockerfile()

Checks the Docker build file.

Checks that:
: \* a name is given and is consistent with the pipeline name

- dependency versions are pinned
- dependency versions are the latest available

#### check_conda_env_yaml()

Checks that the conda environment file is valid.

Checks that:
: \* a name is given and is consistent with the pipeline name

- check that dependency versions are pinned
- dependency versions are the latest available

#### check_cookiecutter_strings()

Look for the string ‘cookiecutter’ in all pipeline files.
Finding it probably means that there has been a copy+paste error from the template.

#### check_docker()

Checks that Dockerfile contains the string `FROM`.

#### check_files_exist()

Checks a given pipeline directory for required files.

Iterates through the pipeline’s directory content and checkmarks files
for presence.
Files that **must** be present:

```default
'nextflow.config',
'nextflow_schema.json',
'Dockerfile',
['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling
'README.md',
'CHANGELOG.md',
'docs/README.md',
'docs/output.md',
'docs/usage.md',
'.github/workflows/branch.yml',
'.github/workflows/ci.yml',
'.github/workflows/linting.yml'
```

Files that _should_ be present:

```default
'main.nf',
'environment.yml',
'conf/base.config',
'.github/workflows/awstest.yml',
'.github/workflows/awsfulltest.yml'
```

Files that _must not_ be present:

```default
'Singularity'
```

Files that _should not_ be present:

```default
'.travis.yml'
```

- **Raises:**
  **An AssertionError if neither nextflow.config** **or** **main.nf found.** –

#### check_licence()

Checks licence file is MIT.

Currently the checkpoints are:
: \* licence file must be long enough (4 or more lines)

- licence contains the string _without restriction_
- licence doesn’t have any placeholder variables

#### check_nextflow_config()

Checks a given pipeline for required config variables.

At least one string in each list must be present for fail and warn.
Any config in config_fail_ifdefined results in a failure.

Uses `nextflow config -flat` to parse pipeline `nextflow.config`
and print all config variables.
NB: Does NOT parse contents of main.nf / nextflow script

#### check_pip_package(dep)

Query PyPi package information.

Sends a HTTP GET request to the PyPi remote API.

- **Parameters:**
  **dep** (_str_) – A PyPi package name.
- **Raises:**
  **A ValueError\*\***,\*\* **if the package name can not be resolved** **or** **the connection timed out.** –

#### check_pipeline_name()

Check whether pipeline name adheres to lower case/no hyphen naming convention

#### check_pipeline_todos()

Go through all template files looking for the string ‘TODO nf-core:’

#### check_readme()

Checks the repository README file for errors.

Currently just checks the badges at the top of the README.

#### check_schema_lint()

Lint the pipeline schema

#### check_schema_params()

Check that the schema describes all flat params in the pipeline

#### check_version_consistency()

Checks container tags versions.

Runs on `process.container` (if set) and `$GITHUB_REF` (if a GitHub Actions release).

Checks that:
: \* the container has a tag

- the version numbers are numeric
- the version numbers are the same as one-another

#### get_results_md()

Function to create a markdown file suitable for posting in a GitHub comment

#### github_comment()

If we are running in a GitHub PR, try to post results as a comment

#### lint_pipeline(release_mode=False)

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
  **If a critical problem is found\*\***,\*\* **an AssertionError is raised.** –

#### save_json_results(json_fn)

Function to dump lint results to a JSON file for downstream use

### nf_core.lint.run_linting(pipeline_dir, release_mode=False, show_passed=False, md_fn=None, json_fn=None)

Runs all nf-core linting checks on a given Nextflow pipeline project
in either release mode or normal mode (default). Returns an object
of type [`PipelineLint`](#nf_core.lint.PipelineLint) after finished.

- **Parameters:**
  - **pipeline_dir** (_str_) – The path to the Nextflow pipeline root directory
  - **release_mode** (_bool_) – Set this to True, if the linting should be run in the release mode.
    See [`PipelineLint`](#nf_core.lint.PipelineLint) for more information.
- **Returns:**
  An object of type [`PipelineLint`](#nf_core.lint.PipelineLint) that contains all the linting results.
