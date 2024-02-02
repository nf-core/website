# nf_core.download

Downloads a nf-core pipeline to the local file system.

### _class_ nf_core.download.DownloadWorkflow(pipeline, release=None, singularity=False, outdir=None)

Downloads a nf-core workflow from Github to the local file system.

Can also download its Singularity container image if required.

- **Parameters:**
  - **pipeline** (_str_) – A nf-core pipeline name.
  - **release** (_str_) – The workflow release version to download, like 1.0. Defaults to None.
  - **singularity** (_bool_) – Flag, if the Singularity container should be downloaded as well. Defaults to False.
  - **outdir** (_str_) – Path to the local download directory. Defaults to None.

#### download_wf_files()

Downloads workflow files from Github to the `self.outdir`.

#### download_workflow()

Starts a nf-core workflow download.

#### fetch_workflow_details(wfs)

Fetches details of a nf-core workflow to download.

- **Parameters:**
  **wfs** ([_nf_core.list.Workflows_](list.md#nf_core.list.Workflows)) – A nf_core.list.Workflows object
- **Raises:**
  **LockupError\*\***,\*\* **if the pipeline can not be found.** –

#### find_container_images()

Find container image names for workflow

#### pull_singularity_image(container)

Uses a local installation of singularity to pull an image from Docker Hub.

- **Parameters:**
  **container** (_str_) – A pipeline’s container name. Usually it is of similar format
  to nfcore/name:dev.
- **Raises:**
  **Various exceptions possible from subprocess execution** **of** **Singularity.** –

#### validate_md5(fname, expected)

Calculates the md5sum for a file on the disk and validate with expected.

- **Parameters:**
  - **fname** (_str_) – Path to a local file.
  - **expected** (_str_) – The expected md5sum.
- **Raises:**
  **IOError\*\***,\*\* **if the md5sum does not match the remote sum.** –
