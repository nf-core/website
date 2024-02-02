# nf_core.sync

Synchronise a pipeline TEMPLATE branch with the template.

### _class_ nf_core.sync.PipelineSync(pipeline_dir, from_branch=None, make_pr=False, gh_repo=None, gh_username=None)

Bases: `object`

Object to hold syncing information and results.

- **Parameters:**
  - **pipeline_dir** (_str_) – The path to the Nextflow pipeline root directory
  - **from_branch** (_str_) – The branch to use to fetch config vars. If not set, will use current active branch
  - **make_pr** (_bool_) – Set this to True to create a GitHub pull-request with the changes
  - **gh_username** (_str_) – GitHub username
  - **gh_repo** (_str_) – GitHub repository name

#### pipeline_dir

Path to target pipeline directory

- **Type:**
  str

#### from_branch

Repo branch to use when collecting workflow variables. Default: active branch.

- **Type:**
  str

#### original_branch

Repo branch that was checked out before we started.

- **Type:**
  str

#### made_changes

Whether making the new template pipeline introduced any changes

- **Type:**
  bool

#### make_pr

Whether to try to automatically make a PR on GitHub.com

- **Type:**
  bool

#### required_config_vars

List of nextflow variables required to make template pipeline

- **Type:**
  list

#### gh_username

GitHub username

- **Type:**
  str

#### gh_repo

GitHub repository name

- **Type:**
  str

#### checkout_template_branch()

Try to check out the origin/TEMPLATE in a new TEMPLATE branch.
If this fails, try to check out an existing local TEMPLATE branch.

#### close_open_pr(pr)

Given a PR API response, add a comment and close.

#### close_open_template_merge_prs()

Get all template merging branches (starting with ‘nf-core-template-merge-‘)
and check for any open PRs from these branches to the self.from_branch
If open PRs are found, add a comment and close them

#### commit_template_changes()

If we have any changes with the new template files, make a git commit

#### create_merge_base_branch()

Create a new branch from the updated TEMPLATE branch
This branch will then be used to create the PR

#### delete_template_branch_files()

Delete all files in the TEMPLATE branch

#### get_wf_config()

Check out the target branch if requested and fetch the nextflow config.
Check that we have the required config variables.

#### inspect_sync_dir()

Takes a look at the target directory for syncing. Checks that it’s a git repo
and makes sure that there are no uncommitted changes.

#### make_pull_request()

Create a pull request to a base branch (default: dev),
from a head branch (default: TEMPLATE)

Returns: An instance of class requests.Response

#### make_template_pipeline()

Delete all files and make a fresh template using the workflow variables

#### push_merge_branch()

Push the newly created merge branch to the remote repository

#### push_template_branch()

If we made any changes, push the TEMPLATE branch to the default remote
and try to make a PR. If we don’t have the auth token, try to figure out a URL
for the PR and print this to the console.

#### reset_target_dir()

Reset the target pipeline directory. Check out the original branch.

#### sync()

Find workflow attributes, create a new template pipeline on TEMPLATE

### _exception_ nf_core.sync.PullRequestException

Bases: `Exception`

Exception raised when there was an error creating a Pull-Request on GitHub.com

### _exception_ nf_core.sync.SyncException

Bases: `Exception`

Exception raised when there was an error with TEMPLATE branch synchronisation
