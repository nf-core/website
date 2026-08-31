# meta_yml

#### `ModuleLint.meta_yml(module: NFCoreComponent, allow_missing: bool = False) → None{:python}`

Lint a `meta.yml` file

Checks that the module has a `meta.yml` file, that it is valid according
to the nf-core JSON schema, and that its contents are consistent with
`main.nf`.

The following checks are performed:

- `meta_yml_exists`: The `meta.yml` file must exist.
- `meta_yml_valid`: The `meta.yml` must be valid according to the JSON
  schema defined in `modules/meta-schema.json` in the nf-core/modules
  repository.
- `meta_name`: The `name` field in `meta.yml` must match (case-insensitive)
  the process name declared in `main.nf`.
- `meta_input`: If `main.nf` declares inputs, they must be listed under
  the `input:` key in `meta.yml`.
- `correct_meta_inputs`: The inputs listed in `meta.yml` must exactly
  match those parsed from `main.nf`. Run `nf-core modules lint --fix`
  to auto-correct.
- `meta_output`: If `main.nf` declares outputs, they must be listed under
  the `output:` key in `meta.yml`.
- `correct_meta_outputs`: The outputs listed in `meta.yml` must exactly
  match those parsed from `main.nf`. Run `nf-core modules lint --fix`
  to auto-correct.
- `has_meta_topics`: If `main.nf` declares topics, `meta.yml` must
  also contain a non-empty `topics:` block. Run
  `nf-core modules lint --fix` to auto-correct.
- `correct_meta_topics`: The topics listed in `meta.yml` must exactly
  match those parsed from `main.nf`. Run `nf-core modules lint --fix`
  to auto-correct.

If the module has inputs or outputs, they are expected to be formatted as:

```groovy
tuple val(foo) path(bar)
val foo
path foo
```

or permutations of the above.
