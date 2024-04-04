# nf\_core.licences

Lists software licences for a given workflow.

### *`class{:python}`*`nf_core.licences.WorkflowLicences(pipeline){:python}`

Bases: `object`

A nf-core workflow licenses collection.

Tries to retrieve the license information from all dependencies
of a given nf-core pipeline.

A condensed overview with license per dependency can be printed out.

* **Parameters:**
  **pipeline** (*str*) – An existing nf-core pipeline name, like nf-core/hlatyping
  or short hlatyping.

#### `fetch_conda_licences(){:python}`

Fetch package licences from Anaconda and PyPi.

#### `get_environment_file(){:python}`

Get the conda environment file for the pipeline

#### `print_licences(){:python}`

Prints the fetched license information.

* **Parameters:**
  **as\_json** (*boolean*) – Prints the information in JSON. Defaults to False.

#### `run_licences(){:python}`

Run the nf-core licences action
