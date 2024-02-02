# nf_core.licences

Lists software licences for a given workflow.

### _class_ nf_core.licences.WorkflowLicences(pipeline)

Bases: `object`

A nf-core workflow licenses collection.

Tries to retrieve the license information from all dependencies
of a given nf-core pipeline.

A condensed overview with license per dependency can be printed out.

- **Parameters:**
  **pipeline** (_str_) – An existing nf-core pipeline name, like nf-core/hlatyping
  or short hlatyping.

#### fetch_conda_licences()

Fetch package licences from Anaconda and PyPi.

#### get_environment_file()

Get the conda environment file for the pipeline

#### print_licences()

Prints the fetched license information.

- **Parameters:**
  **as_json** (_boolean_) – Prints the information in JSON. Defaults to False.

#### run_licences()

Run the nf-core licences action
