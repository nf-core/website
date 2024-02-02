# nf_core.create

Creates a nf-core pipeline matching the current
organization’s specification based on a template.

### _class_ nf_core.create.PipelineCreate(name, description, author, version='1.0dev', no_git=False, force=False, outdir=None)

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

#### download_pipeline_logo(url, img_fn)

Attempt to download a logo from the website. Retry if it fails.

#### git_init_pipeline()

Initialises the new pipeline as a Git repository and submits first commit.

#### init_pipeline()

Creates the nf-core pipeline.

#### make_pipeline_logo()

Fetch a logo for the new pipeline from the nf-core website

#### render_template()

Runs Jinja to create a new nf-core pipeline.
