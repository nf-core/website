# nextflow_config

#### PipelineLint.nextflow_config()

Checks the pipeline configuration for required variables.

All nf-core pipelines are required to be configured with a minimal set of variable
names. This test fails or throws warnings if required variables are not set.

#### NOTE

These config variables must be set in `nextflow.config` or another config
file imported from there. Any variables set in nextflow script files (eg. `main.nf`)
are not checked and will be assumed to be missing.

**The following variables fail the test if missing:**

- `params.outdir`: A directory in which all pipeline results should be saved
- `manifest.name`: The pipeline name. Should begin with `nf-core/`
- `manifest.description`: A description of the pipeline
- `manifest.version`
  - The version of this pipeline. This should correspond to a [GitHub release](https://help.github.com/articles/creating-releases/).
  - If `--release` is set when running `nf-core lint`, the version number must not contain the string `dev`
  - If `--release` is \_not\_ set, the version should end in `dev` (warning triggered if not)
- `manifest.nextflowVersion`
  - The minimum version of Nextflow required to run the pipeline.
  - Should be `>=` or `!>=` and a version number, eg. `manifest.nextflowVersion = '>=0.31.0'` (see [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#scope-manifest))
  - `>=` warns about old versions but tries to run anyway, `!>=` fails for old versions. Only use the latter if you _know_ that the pipeline will certainly fail before this version.
  - This should correspond to the `NXF_VER` version tested by GitHub Actions.
- `manifest.homePage`
  - The homepage for the pipeline. Should be the nf-core GitHub repository URL,
    so beginning with `https://github.com/nf-core/`
- `timeline.enabled`, `trace.enabled`, `report.enabled`, `dag.enabled`
  - The nextflow timeline, trace, report and DAG should be enabled by default (set to `true`)
- `process.cpus`, `process.memory`, `process.time`
  - Default CPUs, memory and time limits for tasks
- `params.input`
  - Input parameter to specify input data, specify this to avoid a warning
  - Typical usage:
    - `params.input`: Input data that is not NGS sequencing data
- `params.custom_config_version`
  > - Should always be set to default value `master`
- `params.custom_config_base`

  > - Should always be set to default value:

  > `https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}`

- `params.validationShowHiddenParams`
  > - Determines whether boilerplate params are showed by schema. Set to `false` by default
- `params.validationSchemaIgnoreParams`
  > - A comma separated string of inputs the schema validation should ignore.

**The following variables throw warnings if missing:**

- `manifest.mainScript`: The filename of the main pipeline script (should be `main.nf`)
- `timeline.file`, `trace.file`, `report.file`, `dag.file`
  - Default filenames for the timeline, trace and report
  - The DAG file path should end with `.svg` (If Graphviz is not installed, Nextflow will generate a `.dot` file instead)

**The following variables are depreciated and fail the test if they are still present:**

- `params.version`: The old method for specifying the pipeline version. Replaced by `manifest.version`
- `params.nf_required_version`: The old method for specifying the minimum Nextflow version. Replaced by `manifest.nextflowVersion`
- `params.container`: The old method for specifying the dockerhub container address. Replaced by `process.container`
- `igenomesIgnore`: Changed to `igenomes_ignore`

**The following Nextflow syntax is depreciated and fails the test if present:**

- Process-level configuration syntax still using the old Nextflow syntax, for example: `process.$fastqc` instead of `process withName:'fastqc'`.

**The configuration should contain the following or the test will fail:**

- A `test` configuration profile should exist.

**The default values in \`\`nextflow.config\`\` should match the default values defined in the \`\`nextflow_schema.json\`\`.**
