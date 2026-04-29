# meta_yml

#### `SubworkflowLint.meta_yml(subworkflow, allow_missing: bool = False){:python}`

Lint a `meta.yml` file

The lint test checks that the subworkflow has
a `meta.yml` file and that it follows the
JSON schema defined in the `subworkflows/yaml-schema.json`
file in the nf-core/modules repository.

In addition it checks that the subworkflow name
and subworkflow input is consistent between the
`meta.yml` and the `main.nf`.

Checks that all input and output channels are specified in `meta.yml`.
Checks that all included components in `main.nf` are specified in `meta.yml`.

The following checks are performed:

- `meta_yml_exists`: The `meta.yml` file must exist.
- `meta_yml_valid`: The `meta.yml` must be valid according to the JSON
  schema defined at <https://raw.githubusercontent.com/nf-core/subworkflows/master/modules/environment-schema.json>.
- `meta_input`: All input channels declared in `main.nf` must be listed
  under the `input:` key in `meta.yml`.
- `meta_output`: All output channels declared in `main.nf` must be listed
  under the `output:` key in `meta.yml`.
- `meta_name`: The `name` field in `meta.yml` must match the workflow
  name declared in `main.nf`.
- `meta_include`: All modules and subworkflows included in `main.nf` must
  be listed under the `components:` key in `meta.yml`.
- `meta_modules_deprecated`: The deprecated `modules:` section must not
  be present in `meta.yml`; use `components:` instead.
