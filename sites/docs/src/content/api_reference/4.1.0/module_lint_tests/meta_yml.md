# meta\_yml

#### `ModuleLint.meta_yml(module: NFCoreComponent, allow_missing: bool = False) → None{:python}`

Lint a `meta.yml` file

Checks that the module has a `meta.yml` file, that it is valid according
to the nf-core JSON schema, and that its contents are consistent with
`main.nf`.

The following checks are performed:

## meta\_yml\_exists

The `meta.yml` file must exist.

## meta\_yml\_valid

The `meta.yml` must be valid according to the JSON
schema defined in `modules/meta-schema.json` in the nf-core/modules
repository.

## meta\_name

The `name` field in `meta.yml` must match (case-insensitive)
the process name declared in `main.nf`.

## meta\_input

If `main.nf` declares inputs, they must be listed under
the `input:` key in `meta.yml`.

## correct\_meta\_inputs

The inputs listed in `meta.yml` must exactly
match those parsed from `main.nf`. Run `nf-core modules lint --fix`
to auto-correct.

## meta\_output

If `main.nf` declares outputs, they must be listed under
the `output:` key in `meta.yml`.

## correct\_meta\_outputs

The outputs listed in `meta.yml` must exactly
match those parsed from `main.nf`. Run `nf-core modules lint --fix`
to auto-correct.

## has\_meta\_topics

If `main.nf` declares topics, `meta.yml` must
also contain a non-empty `topics:` block. Run
`nf-core modules lint --fix` to auto-correct.

## correct\_meta\_topics

The topics listed in `meta.yml` must exactly
match those parsed from `main.nf`. Run `nf-core modules lint --fix`
to auto-correct.

## has\_meta\_containers

If `main.nf` declares containers, `meta.yml`
must also contain a non-empty `containers:` block. Run
`nf-core modules lint --fix` to auto-correct.

## correct\_meta\_containers

The containers listed in `meta.yml` must
exactly match those parsed from `main.nf`. Run
`nf-core modules lint --fix` to auto-correct.

If the module has inputs or outputs, they are expected to be formatted as:

```groovy
tuple val(foo) path(bar)
val foo
path foo
```

or permutations of the above.
