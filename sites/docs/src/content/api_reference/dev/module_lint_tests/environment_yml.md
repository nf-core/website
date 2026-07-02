# environment\_yml

#### `ModuleLint.environment_yml(module: NFCoreComponent, allow_missing: bool = False, fix_version: bool = False, progress_bar: Progress | None = None) → None{:python}`

Lint an `environment.yml` file.

The lint test checks that the `dependencies` section
in the environment.yml file is valid YAML and that it
is sorted alphabetically.

The following checks are performed:

## environment\_yml\_exists

The `environment.yml` file must exist, unless the
module has a Dockerfile. If neither `environment.yml` nor a Dockerfile is
present, the test fails.

## environment\_yml\_dockerfile\_conflict

A module must not have both an
`environment.yml` and a Dockerfile; Remove the Dockerfile if present.

## environment\_yml\_valid

The `environment.yml` must be valid according to
the JSON schema defined at <https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json>.

## environment\_yml\_sorted

The dependencies listed in `environment.yml`
must be sorted alphabetically. If they are not, they will be sorted
automatically.

## bioconda\_version

The version specified for each bioconda package must
be a valid, known version available in the bioconda channel.

## bioconda\_latest

Each bioconda package should be pinned to the latest
available version. A warning is raised (or the version updated if
`fix_version` is set) when a newer version exists.
