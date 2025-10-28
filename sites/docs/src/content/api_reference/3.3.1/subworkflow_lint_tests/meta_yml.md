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
