# meta_yml

#### ModuleLint.meta_yml(module)

Lint a `meta.yml` file

The lint test checks that the module has
a `meta.yml` file and that it follows the
JSON schema defined in the `modules/yaml-schema.json`
file in the nf-core/modules repository.

In addition it checks that the module name
and module input is consistent between the
`meta.yml` and the `main.nf`.
