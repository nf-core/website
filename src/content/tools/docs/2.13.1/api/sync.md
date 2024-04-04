# nf\_core.sync

Synchronise a pipeline TEMPLATE branch with the template.

### *`class{:python}`*`nf_core.sync.PipelineSync(pipeline_dir, from_branch=None, make_pr=False, gh_repo=None, gh_username=None, template_yaml_path=None){:python}`

Bases: `object`

Object to hold syncing information and results.

* **Parameters:**
  * **pipeline\_dir** (*str*) – The path to the Nextflow pipeline root directory
  * **from\_branch** (*str*) – The branch to use to fetch config vars. If not set, will use current active branch
  * **make\_pr** (*bool*) – Set this to True to create a GitHub pull-request with the changes
  * **gh\_username** (*str*) – GitHub username
  * **gh\_repo** (*str*) – GitHub repository name
  * **template\_yaml\_path** (*str*) – Path to template.yml file for pipeline creation settings. DEPRECATED

#### `pipeline_dir{:python}`

Path to target pipeline directory

* **Type:**
  str

#### `from_branch{:python}`

Repo branch to use when collecting workflow variables. Default: active branch.

* **Type:**
  str

#### `original_branch{:python}`

Repo branch that was checked out before we started.

* **Type:**
  str

#### `made_changes{:python}`

Whether making the new template pipeline introduced any changes

* **Type:**
  bool

#### `make_pr{:python}`

Whether to try to automatically make a PR on GitHub.com

* **Type:**
  bool

#### `required_config_vars{:python}`

List of nextflow variables required to make template pipeline

* **Type:**
  list

#### `gh_username{:python}`

GitHub username

* **Type:**
  str

#### `gh_repo{:python}`

GitHub repository name

* **Type:**
  str

#### `checkout_template_branch(){:python}`

Try to check out the origin/TEMPLATE in a new TEMPLATE branch.
If this fails, try to check out an existing local TEMPLATE branch.

#### `close_open_pr(pr){:python}`

Given a PR API response, add a comment and close.

#### `close_open_template_merge_prs(){:python}`

Get all template merging branches (starting with ‘nf-core-template-merge-‘)
and check for any open PRs from these branches to the self.from\_branch
If open PRs are found, add a comment and close them

#### `commit_template_changes(){:python}`

If we have any changes with the new template files, make a git commit

#### `create_merge_base_branch(){:python}`

Create a new branch from the updated TEMPLATE branch
This branch will then be used to create the PR

#### `delete_template_branch_files(){:python}`

Delete all files in the TEMPLATE branch

#### `get_wf_config(){:python}`

Check out the target branch if requested and fetch the nextflow config.
Check that we have the required config variables.

#### `inspect_sync_dir(){:python}`

Takes a look at the target directory for syncing. Checks that it’s a git repo
and makes sure that there are no uncommitted changes.

#### `make_pull_request(){:python}`

Create a pull request to a base branch (default: dev),
from a head branch (default: TEMPLATE)

Returns: An instance of class requests.Response

#### `make_template_pipeline(){:python}`

Delete all files and make a fresh template using the workflow variables

#### `push_merge_branch(){:python}`

Push the newly created merge branch to the remote repository

#### `push_template_branch(){:python}`

If we made any changes, push the TEMPLATE branch to the default remote
and try to make a PR. If we don’t have the auth token, try to figure out a URL
for the PR and print this to the console.

#### `reset_target_dir(){:python}`

Reset the target pipeline directory. Check out the original branch.

#### `sync(){:python}`

Find workflow attributes, create a new template pipeline on TEMPLATE

### *`exception{:python}`*`nf_core.sync.PullRequestExceptionError{:python}`

Bases: `Exception`

Exception raised when there was an error creating a Pull-Request on GitHub.com

### *`exception{:python}`*`nf_core.sync.SyncExceptionError{:python}`

Bases: `Exception`

Exception raised when there was an error with TEMPLATE branch synchronisation
