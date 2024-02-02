# nf_core.create

Creates a nf-core pipeline matching the current
organization’s specification based on a template.

### _class_ nf_core.create.PipelineCreate(name, description, author, version='1.0dev', no_git=False, force=False, outdir=None, template_yaml_path=None, plain=False, default_branch=None)

Bases: `object`

Creates a nf-core pipeline a la carte from the nf-core best-practice template.

- **Parameters:**
  - **name** (_str_) – Name for the pipeline.
  - **description** (_str_) – Description for the pipeline.
  - **author** (_str_) – Authors name of the pipeline.
  - **version** (_str_) – Version flag. Semantic versioning only. Defaults to 1.0dev.
  - **no_git** (_bool_) – Prevents the creation of a local Git repository for the pipeline. Defaults to False.
  - **force** (_bool_) – Overwrites a given workflow directory with the same name. Defaults to False.
    May the force be with you.
  - **outdir** (_str_) – Path to the local output directory.
  - **template_yaml_path** (_str_) – Path to template.yml file for pipeline creation settings.
  - **plain** (_bool_) – If true the Git repository will be initialized plain.
  - **default_branch** (_str_) – Specifies the –initial-branch name.

#### create_param_dict(name, description, author, version, template_yaml_path, plain)

Creates a dictionary of parameters for the new pipeline.

- **Parameters:**
  **template_yaml_path** (_str_) – Path to YAML file containing template parameters.

#### customize_template(template_areas)

Customizes the template parameters.

- **Parameters:**
  - **name** (_str_) – Name for the pipeline.
  - **description** (_str_) – Description for the pipeline.
  - **author** (_str_) – Authors name of the pipeline.

#### download_pipeline_logo(url, img_fn)

Attempt to download a logo from the website. Retry if it fails.

#### fix_linting()

Updates the .nf-core.yml with linting configurations
for a customized pipeline.

#### get_param(param_name, passed_value, template_yaml, template_yaml_path)

#### git_init_pipeline()

Initialises the new pipeline as a Git repository and submits first commit.

- **Raises:**
  **UserWarning** – if Git default branch is set to ‘dev’ or ‘TEMPLATE’.

#### init_pipeline()

Creates the nf-core pipeline.

#### make_pipeline_logo()

Fetch a logo for the new pipeline from the nf-core website

#### prompt_wf_author()

#### prompt_wf_description()

#### prompt_wf_name()

#### remove_nf_core_in_bug_report_template()

Remove the field mentioning nf-core documentation
in the github bug report template

#### render_template()

Runs Jinja to create a new nf-core pipeline.

#### update_nextflow_schema()

Removes unused parameters from the nextflow schema.
