# environment_yml

#### `ModuleLint.environment_yml(module: NFCoreComponent, allow_missing: bool = False) → None{:python}`

Lint an `environment.yml` file.

The lint test checks that the `dependencies` section
in the environment.yml file is valid YAML and that it
is sorted alphabetically.

The following checks are performed:

- `environment_yml_exists`: The `environment.yml` file must exist if it is
  referenced in `main.nf`.
- `environment_yml_valid`: The `environment.yml` must be valid according to
  the JSON schema defined at <https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json>.
- `environment_yml_sorted`: The dependencies listed in `environment.yml`
  must be sorted alphabetically. If they are not, they will be sorted
  automatically.
