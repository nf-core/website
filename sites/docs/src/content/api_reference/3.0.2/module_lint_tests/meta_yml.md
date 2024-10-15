# meta_yml

#### `ModuleLint.meta_yml(module: NFCoreComponent) → None{:python}`

Lint a `meta.yml` file

The lint test checks that the module has
a `meta.yml` file and that it follows the
JSON schema defined in the `modules/meta-schema.json`
file in the nf-core/modules repository.

In addition it checks that the module name
and module input is consistent between the
`meta.yml` and the `main.nf`.

If the module has inputs or outputs, they are expected to be
formatted as:

```groovy
tuple val(foo) path(bar)
val foo
path foo
```

or permutations of the above.

- **Parameters:**
  - **module_lint_object** (_ComponentLint_) – The lint object for the module
  - **module** (_NFCoreComponent_) – The module to lint
