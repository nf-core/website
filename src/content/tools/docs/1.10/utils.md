# nf_core.utils

Common utility functions for the nf-core python package.

### nf_core.utils.check_if_outdated(current_version=None, remote_version=None, source_url='https://nf-co.re/tools_version')

Check if the current version of nf-core is outdated

### nf_core.utils.fetch_wf_config(wf_path)

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

- **Parameters:**
  **wf_path** (_str_) – Nextflow workflow file system path.
- **Returns:**
  Workflow configuration settings.
- **Return type:**
  dict

### nf_core.utils.poll_nfcore_web_api(api_url, post_data=None)

Poll the nf-core website API

Takes argument api_url for URL

Expects API reponse to be valid JSON and contain a top-level ‘status’ key.

### nf_core.utils.setup_requests_cachedir()

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .nfcore_cache subdir.

### nf_core.utils.wait_cli_function(poll_func, poll_every=20)

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

- **Parameters:**
  - **poll_func** (_function_) – Function to call
  - **poll_every** (_int_) – How many tenths of a second to wait between function calls. Default: 20.
- **Returns:**
  None. Just sits in an infite loop until the function returns True.
