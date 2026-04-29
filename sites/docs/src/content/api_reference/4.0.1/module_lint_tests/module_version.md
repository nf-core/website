# module_version

#### `ModuleLint.module_version(module: NFCoreComponent){:python}`

Verifies that the module has a version specified in the `modules.json` file

It checks whether the module has an entry in the `modules.json` file
containing a commit SHA. If that is true, it verifies that there are no
newer version of the module available.

The following checks are performed:

- `git_sha`: The module must have a `git_sha` entry in `modules.json`.
- `module_version`: The module version must match the latest commit in the
  remote repository. A warning is issued if a newer version is available.
