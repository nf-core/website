# nf\_core.create

Creates a nf-core pipeline matching the current
organization’s specification based on a template.

### *`class{:python}`*`nf_core.create.PipelineCreate(name, description, author, version='1.0dev', no_git=False, force=False, outdir=None, template_yaml_path=None, plain=False, default_branch=None){:python}`

Bases: `object`

Creates a nf-core pipeline a la carte from the nf-core best-practice template.

* **Parameters:**
  * **name** (*str*) – Name for the pipeline.
  * **description** (*str*) – Description for the pipeline.
  * **author** (*str*) – Authors name of the pipeline.
  * **version** (*str*) – Version flag. Semantic versioning only. Defaults to 1.0dev.
  * **no\_git** (*bool*) – Prevents the creation of a local Git repository for the pipeline. Defaults to False.
  * **force** (*bool*) – Overwrites a given workflow directory with the same name. Defaults to False.
    May the force be with you.
  * **outdir** (*str*) – Path to the local output directory.
  * **template\_yaml\_path** (*str*) – Path to template.yml file for pipeline creation settings.
  * **plain** (*bool*) – If true the Git repository will be initialized plain.
  * **default\_branch** (*str*) – Specifies the –initial-branch name.

#### `create_param_dict(name, description, author, version, template_yaml_path, plain, pipeline_dir){:python}`

Creates a dictionary of parameters for the new pipeline.

* **Parameters:**
  * **name** (*str*) – Name for the pipeline.
  * **description** (*str*) – Description for the pipeline.
  * **author** (*str*) – Authors name of the pipeline.
  * **version** (*str*) – Version flag.
  * **template\_yaml\_path** (*str*) – Path to YAML file containing template parameters.
  * **plain** (*bool*) – If true the pipeline template will be initialized plain, without customisation.
  * **pipeline\_dir** (*str*) – Path to the pipeline directory.

#### `customize_template(template_areas){:python}`

Customizes the template parameters.

* **Parameters:**
  **template\_areas** (*list<str>*) – List of available template areas to skip.

#### `fix_linting(){:python}`

Updates the .nf-core.yml with linting configurations
for a customized pipeline.

#### `get_param(param_name, passed_value, template_yaml, template_yaml_path){:python}`

#### `git_init_pipeline(){:python}`

Initialises the new pipeline as a Git repository and submits first commit.

* **Raises:**
  **UserWarning** – if Git default branch is set to ‘dev’ or ‘TEMPLATE’.

#### `init_pipeline(){:python}`

Creates the nf-core pipeline.

#### `make_pipeline_logo(){:python}`

Fetch a logo for the new pipeline from the nf-core website

#### `prompt_wf_author(){:python}`

#### `prompt_wf_description(){:python}`

#### `prompt_wf_name(){:python}`

#### `remove_nf_core_in_bug_report_template(){:python}`

Remove the field mentioning nf-core documentation
in the github bug report template

#### `render_template(){:python}`

Runs Jinja to create a new nf-core pipeline.

#### `update_nextflow_schema(){:python}`

Removes unused parameters from the nextflow schema.
