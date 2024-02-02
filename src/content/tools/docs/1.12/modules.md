# nf_core.modules

Code to handle DSL2 module imports from a GitHub repository

### _class_ nf_core.modules.ModulesRepo(repo='nf-core/modules', branch='master')

Bases: `object`

An object to store details about the repository being used for modules.

Used by the nf-core modules top-level command with -r and -b flags,
so that this can be used in the same way by all sucommands.

### _class_ nf_core.modules.PipelineModules

Bases: `object`

#### check_modules()

#### download_gh_file(dl_filename, api_url)

Download a file from GitHub using the GitHub API

- **Parameters:**
  - **dl_filename** (_string_) – Path to save file to
  - **api_url** (_string_) – GitHub API URL for file
- **Raises:**
  **If a problem\*\***,\*\* **raises an error** –

#### get_module_file_urls(module)

Fetch list of URLs for a specific module

Takes the name of a module and iterates over the GitHub repo file tree.
Loops over items that are prefixed with the path ‘software/<module_name>’ and ignores
anything that’s not a blob. Also ignores the test/ subfolder.

Returns a dictionary with keys as filenames and values as GitHub API URIs.
These can be used to then download file contents.

- **Parameters:**
  **module** (_string_) – Name of module for which to fetch a set of URLs
- **Returns:**
  Set of files and associated URLs as follows:
  {
  : ‘software/fastqc/main.nf’: ‘[https://api.github.com/repos/nf-core/modules/git/blobs/65ba598119206a2b851b86a9b5880b5476e263c3](https://api.github.com/repos/nf-core/modules/git/blobs/65ba598119206a2b851b86a9b5880b5476e263c3)’,
  ‘software/fastqc/meta.yml’: ‘[https://api.github.com/repos/nf-core/modules/git/blobs/0d5afc23ba44d44a805c35902febc0a382b17651](https://api.github.com/repos/nf-core/modules/git/blobs/0d5afc23ba44d44a805c35902febc0a382b17651)’

  }

- **Return type:**
  dict

#### get_modules_file_tree()

Fetch the file list from the repo, using the GitHub API

Sets self.modules_file_tree
: self.modules_current_hash
self.modules_avail_module_names

#### install(module)

#### list_modules()

Get available module names from GitHub tree for repo
and print as list to stdout

#### remove(module)

#### update(module, force=False)
