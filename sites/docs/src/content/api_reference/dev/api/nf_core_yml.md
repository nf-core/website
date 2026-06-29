# .nf-core.yml configuration

The `.nf-core.yml` file at the root of any nf-core repository controls how nf-core
tools behaves for that repository. It is read by `load_tools_config()` and validated
against the `NFCoreYamlConfig` Pydantic model.

## Minimal examples

Pipeline repository:

```yaml
repository_type: pipeline
nf_core_version: "4.0.2"
```

Modules repository:

```yaml
repository_type: modules
nf_core_version: "4.0.2"
container-registry:
  - community.wave.seqera.io/library/
```

## Schema

### `Top-level{:python}``.nf-core.yml``schema ({:python}``NFCoreYamlConfig``){:python}`

### _`pydantic model{:python}`_`nf_core.utils.NFCoreYamlConfig{:python}`

Bases: `BaseModel`

.nf-core.yml configuration file schema

- **Fields:**
  - [`bump_version (dict[str, bool] | None)`](utils#nf_core.utils.NFCoreYamlConfig.bump_version)
  - [`container_registry (list[str] | None)`](utils#nf_core.utils.NFCoreYamlConfig.container_registry)
  - [`lint (nf_core.utils.NFCoreYamlLintConfig | None)`](utils#nf_core.utils.NFCoreYamlConfig.lint)
  - [`nf_core_version (str | None)`](utils#nf_core.utils.NFCoreYamlConfig.nf_core_version)
  - [`org_path (str | None)`](utils#nf_core.utils.NFCoreYamlConfig.org_path)
  - [`repository_type (Literal['pipeline', 'modules'] | None)`](utils#nf_core.utils.NFCoreYamlConfig.repository_type)
  - [`template (nf_core.utils.NFCoreTemplateConfig | None)`](utils#nf_core.utils.NFCoreYamlConfig.template)
  - [`update (dict[str, str | bool | dict[str, str | dict[str, str | bool]]] | None)`](utils#nf_core.utils.NFCoreYamlConfig.update)

#### _`field{:python}`_`bump_version{:python}`_: dict\[str, bool] | None_`{:python}`_= None_

Disable bumping of the version for a module/subworkflow (when repository_type is modules). See <https://nf-co.re/docs/nf-core-tools/modules/bump-versions> for more information.

#### _`field{:python}`_`container_registry{:python}`_: list\[str] | None_`{:python}`_= None_`{:python}`_(alias 'container-registry')_

Additional container registry prefixes allowed when linting container directives.

#### _`field{:python}`_`lint{:python}`_: [NFCoreYamlLintConfig](utils#nf_core.utils.NFCoreYamlLintConfig) | None_`{:python}`_= None_

Pipeline linting configuration, see <https://nf-co.re/docs/nf-core-tools/pipelines/lint#linting-config> for examples and documentation

#### _`field{:python}`_`nf_core_version{:python}`_: str | None_`{:python}`_= None_

Version of nf-core/tools used to create/update the pipeline

#### _`field{:python}`_`org_path{:python}`_: str | None_`{:python}`_= None_

Path to the organisation’s modules repository (used for modules repo_type only)

#### _`field{:python}`_`repository_type{:python}`_: Literal\['pipeline', 'modules'] | None_`{:python}`_= None_

Type of repository

#### _`field{:python}`_`template{:python}`_: [NFCoreTemplateConfig](utils#nf_core.utils.NFCoreTemplateConfig) | None_`{:python}`_= None_

Pipeline template configuration

#### _`field{:python}`_`update{:python}`_: dict\[str, str | bool | dict\[str, str | dict\[str, str | bool]]] | None_`{:python}`_= None_

Disable updating specific modules/subworkflows (when repository_type is pipeline). See <https://nf-co.re/docs/nf-core-tools/modules/update> for more information.

### `lint``block ({:python}``NFCoreYamlLintConfig``){:python}`

### _`pydantic model{:python}`_`nf_core.utils.NFCoreYamlLintConfig{:python}`

Bases: `BaseModel`

schema for linting config in .nf-core.yml should cover:

```yaml
files_unchanged:
  - .github/workflows/branch.yml
modules_config: False
modules_config:
  - fastqc
# merge_markers: False
merge_markers:
  - docs/my_pdf.pdf
nextflow_config: False
nextflow_config:
  - manifest.name
  - config_defaults:
      - params.annotation_db
      - params.multiqc_comment_headers
      - params.custom_table_headers
# multiqc_config: False
multiqc_config:
  - report_section_order
  - report_comment
files_exist:
  - CITATIONS.md
template_strings: False
template_strings:
  - docs/my_pdf.pdf
nfcore_components: False
# nf_test_content: False
nf_test_content:
  - tests/<test_name>.nf.test
  - tests/nextflow.config
  - nf-test.config
```

- **Fields:**
  - [`actions_awsfulltest (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.actions_awsfulltest)
  - [`actions_awstest (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.actions_awstest)
  - [`actions_nf_test (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.actions_nf_test)
  - [`actions_schema_validation (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.actions_schema_validation)
  - [`base_config (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.base_config)
  - [`container_configs (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.container_configs)
  - [`files_exist (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.files_exist)
  - [`files_unchanged (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.files_unchanged)
  - [`included_configs (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.included_configs)
  - [`local_component_structure (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.local_component_structure)
  - [`merge_markers (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.merge_markers)
  - [`modules_config (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.modules_config)
  - [`modules_json (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.modules_json)
  - [`modules_structure (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.modules_structure)
  - [`multiqc_config (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.multiqc_config)
  - [`nextflow_config (bool | list[str | dict[str, list[str]]] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.nextflow_config)
  - [`nf_test_content (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.nf_test_content)
  - [`nfcore_components (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.nfcore_components)
  - [`nfcore_yml (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.nfcore_yml)
  - [`pipeline_if_empty_null (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.pipeline_if_empty_null)
  - [`pipeline_name_conventions (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.pipeline_name_conventions)
  - [`pipeline_todos (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.pipeline_todos)
  - [`plugin_includes (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.plugin_includes)
  - [`readme (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.readme)
  - [`rocrate_readme_sync (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.rocrate_readme_sync)
  - [`schema_description (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.schema_description)
  - [`schema_lint (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.schema_lint)
  - [`schema_params (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.schema_params)
  - [`system_exit (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.system_exit)
  - [`template_strings (bool | list[str] | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.template_strings)
  - [`version_consistency (bool | None)`](utils#nf_core.utils.NFCoreYamlLintConfig.version_consistency)

#### _`field{:python}`_`actions_awsfulltest{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to run full tests on AWS

#### _`field{:python}`_`actions_awstest{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to run tests on AWS

#### _`field{:python}`_`actions_nf_test{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to use GitHub Actions CI

#### _`field{:python}`_`actions_schema_validation{:python}`_: bool | None_`{:python}`_= None_

Lint GitHub Action workflow files with schema

#### _`field{:python}`_`base_config{:python}`_: bool | None_`{:python}`_= None_

Lint base.config file

#### _`field{:python}`_`container_configs{:python}`_: bool | None_`{:python}`_= None_

Lint that container configuration files in conf/ are up to date

#### _`field{:python}`_`files_exist{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of files that can not exist

#### _`field{:python}`_`files_unchanged{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of files that should not be changed

#### _`field{:python}`_`included_configs{:python}`_: bool | None_`{:python}`_= None_

Lint for included configs

#### _`field{:python}`_`local_component_structure{:python}`_: bool | None_`{:python}`_= None_

Lint local components use correct structure mirroring remote

#### _`field{:python}`_`merge_markers{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of files that should not contain merge markers

#### _`field{:python}`_`modules_config{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of modules that should not be changed

#### _`field{:python}`_`modules_json{:python}`_: bool | None_`{:python}`_= None_

Lint modules.json file

#### _`field{:python}`_`modules_structure{:python}`_: bool | None_`{:python}`_= None_

Lint modules structure

#### _`field{:python}`_`multiqc_config{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of MultiQC config options that be changed

#### _`field{:python}`_`nextflow_config{:python}`_: bool | list\[str | dict\[str, list\[str]]] | None_`{:python}`_= None_

List of Nextflow config files that should not be changed

#### _`field{:python}`_`nf_test_content{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of nf-test content that should not be changed

#### _`field{:python}`_`nfcore_components{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to use nf-core modules and subworkflows

#### _`field{:python}`_`nfcore_yml{:python}`_: bool | None_`{:python}`_= None_

Lint nf-core.yml

#### _`field{:python}`_`pipeline_if_empty_null{:python}`_: bool | None_`{:python}`_= None_

Lint for ifEmpty(null) statements

#### _`field{:python}`_`pipeline_name_conventions{:python}`_: bool | None_`{:python}`_= None_

Lint for pipeline name conventions

#### _`field{:python}`_`pipeline_todos{:python}`_: bool | None_`{:python}`_= None_

Lint for TODOs statements

#### _`field{:python}`_`plugin_includes{:python}`_: bool | None_`{:python}`_= None_

Lint for nextflow plugin

#### _`field{:python}`_`readme{:python}`_: bool | list\[str] | None_`{:python}`_= None_

Lint the README.md file

#### _`field{:python}`_`rocrate_readme_sync{:python}`_: bool | None_`{:python}`_= None_

Lint for README.md and rocrate.json sync

#### _`field{:python}`_`schema_description{:python}`_: bool | None_`{:python}`_= None_

Check that every parameter in the schema has a description.

#### _`field{:python}`_`schema_lint{:python}`_: bool | None_`{:python}`_= None_

Lint nextflow_schema.json file

#### _`field{:python}`_`schema_params{:python}`_: bool | None_`{:python}`_= None_

Lint schema for all params

#### _`field{:python}`_`system_exit{:python}`_: bool | None_`{:python}`_= None_

Lint for System.exit calls in groovy/nextflow code

#### _`field{:python}`_`template_strings{:python}`_: bool | list\[str] | None_`{:python}`_= None_

List of files that can contain template strings

#### _`field{:python}`_`version_consistency{:python}`_: bool | None_`{:python}`_= None_

Lint for version consistency

### `template``block ({:python}``NFCoreTemplateConfig``){:python}`

### _`pydantic model{:python}`_`nf_core.utils.NFCoreTemplateConfig{:python}`

Bases: `BaseModel`

Template configuration schema

- **Fields:**
  - [`author (str | None)`](utils#nf_core.utils.NFCoreTemplateConfig.author)
  - [`description (str | None)`](utils#nf_core.utils.NFCoreTemplateConfig.description)
  - [`force (bool | None)`](utils#nf_core.utils.NFCoreTemplateConfig.force)
  - [`is_nfcore (bool | None)`](utils#nf_core.utils.NFCoreTemplateConfig.is_nfcore)
  - [`name (str | None)`](utils#nf_core.utils.NFCoreTemplateConfig.name)
  - [`org (str | None)`](utils#nf_core.utils.NFCoreTemplateConfig.org)
  - [`outdir (str | pathlib.Path | None)`](utils#nf_core.utils.NFCoreTemplateConfig.outdir)
  - [`skip_features (list | None)`](utils#nf_core.utils.NFCoreTemplateConfig.skip_features)
  - [`version (str | None)`](utils#nf_core.utils.NFCoreTemplateConfig.version)
- **Validators:**
  - [`outdir_to_str`](utils#nf_core.utils.NFCoreTemplateConfig.outdir_to_str) » [`outdir`](utils#nf_core.utils.NFCoreTemplateConfig.outdir)

#### _`field{:python}`_`author{:python}`_: str | None_`{:python}`_= None_

Pipeline author

#### _`field{:python}`_`description{:python}`_: str | None_`{:python}`_= None_

Pipeline description

#### _`field{:python}`_`force{:python}`_: bool | None_`{:python}`_= True_

Force overwrite of existing files

#### _`field{:python}`_`is_nfcore{:python}`_: bool | None_`{:python}`_= None_

Whether the pipeline is an nf-core pipeline.

#### _`field{:python}`_`name{:python}`_: str | None_`{:python}`_= None_

Pipeline name

#### _`field{:python}`_`org{:python}`_: str | None_`{:python}`_= None_

Organisation name

#### _`field{:python}`_`outdir{:python}`_: str | Path | None_`{:python}`_= None_

Output directory

- **Validated by:**
  - [`outdir_to_str`](utils#nf_core.utils.NFCoreTemplateConfig.outdir_to_str)

#### _`field{:python}`_`skip_features{:python}`_: list | None_`{:python}`_= None_

Skip features. See <https://nf-co.re/docs/nf-core-tools/pipelines/create> for a list of features.

#### _`field{:python}`_`version{:python}`_: str | None_`{:python}`_= None_

Pipeline version
