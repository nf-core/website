# nf_core.list

Lists available nf-core pipelines and versions.

### _class_ nf_core.list.LocalWorkflow(name)

Class to handle local workflows pulled by nextflow

#### get_local_nf_workflow_details()

Get full details about a local cached workflow

### _class_ nf_core.list.RemoteWorkflow(data)

A information container for a remote workflow.

- **Parameters:**
  **data** (_dict_) – workflow information as they are retrieved from the GitHub repository REST API request
  ([https://developer.github.com/v3/repos/#get](https://developer.github.com/v3/repos/#get)).

### _class_ nf_core.list.Workflows(filter_by=None, sort_by='release', show_archived=False)

Workflow container class.

Is used to collect local and remote nf-core pipelines. Pipelines
can be sorted, filtered and compared.

- **Parameters:**
  - **filter_by** (_list_) – A list of strings that can be used for filtering.
  - **sort_by** (_str_) – workflows can be sorted by keywords. Keyword must be one of
    release (default), name, stars.

#### compare_remote_local()

Matches local to remote workflows.

If a matching remote workflow is found, the local workflow’s Git commit hash is compared
with the latest one from remote.

A boolean flag in `RemoteWorkflow.local_is_latest` is set to True, if the local workflow
is the latest.

#### filtered_workflows()

Filters remote workflows for keywords.

- **Returns:**
  Filtered remote workflows.
- **Return type:**
  list

#### get_local_nf_workflows()

Retrieves local Nextflow workflows.

Local workflows are stored in `self.local_workflows` list.

#### get_remote_workflows()

Retrieves remote workflows from [nf-co.re](https://nf-co.re).

Remote workflows are stored in `self.remote_workflows` list.

#### print_json()

Dump JSON of all parsed information

#### print_summary()

Prints a summary of all pipelines.

### nf_core.list.get_local_wf(workflow, revision=None)

Check if this workflow has a local copy and use nextflow to pull it if not

### nf_core.list.list_workflows(filter_by=None, sort_by='release', as_json=False, show_archived=False)

Prints out a list of all nf-core workflows.

- **Parameters:**
  - **filter_by** (_list_) – A list of strings that can be used for filtering.
  - **sort_by** (_str_) – workflows can be sorted by keywords. Keyword must be one of
    release (default), name, stars.
  - **as_json** (_boolean_) – Set to true, if the lists should be printed in JSON.

### nf_core.list.pretty_date(time)

Transforms a datetime object or a int() Epoch timestamp into a
pretty string like ‘an hour ago’, ‘Yesterday’, ‘3 months ago’,
‘just now’, etc

Based on [https://stackoverflow.com/a/1551394/713980](https://stackoverflow.com/a/1551394/713980)
Adapted by sven1103
