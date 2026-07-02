# main\_nf

#### `ModuleLint.main_nf(module: NFCoreComponent, fix_version: bool, registry: tuple[str, ...], progress_bar: Progress) → tuple[list[str], list[str]]{:python}`

Lint a `main.nf` module file

Can also be used to lint local module files,
in which case failures will be reported as warnings.

The following checks are performed:

## main\_nf\_module\_granularity

The module must represent a single command
as `<tool>` or single subcommand with distinct functionality as
`<tool/subtool>`.

## main\_nf\_exists

The `main.nf` file must exist.

## deprecated\_dsl2

The file must not contain deprecated DSL2 identifiers
(`initOptions`, `saveFiles`, `getSoftwareName`, `getProcessName`,
`publishDir`).

## main\_nf\_script\_outputs

The process must have an `output:` block.

## main\_nf\_container

When both a `singularity` and a `docker` container
are specified, their tags must reference the same software version. A warning
is issued if they do not match. Modules using the newer docker-only format
(no singularity container) skip this check.

## singularity\_tag

A Singularity container must be resolvable via
`nextflow inspect -profile singularity`. The check fails if none can be
resolved or if it falls back to a docker container that has an automatic
singularity equivalent, i.e., `quay.io/biocontainers/` or
`community.wave.seqera.io/`.
It is skipped for modules listed under `singularity` in `.github/skip_nf_test.json`.

## oras\_singularity\_tag

The resolved Singularity container must not be
served over the `oras://` scheme; it should be a plain `https://` URL.
The check fails if an `oras://` container is used.

## main\_nf\_script\_shell

Exactly one of `script:`, `shell:`, or `exec:`
blocks must be present.

## main\_nf\_shell\_template

If a `shell:` block is used, it must call
a `template`.

## main\_nf\_meta\_output

If `meta` is present in the module inputs, it
must also appear in at least one output channel.

## main\_nf\_version\_topic

The module should emit software versions using
a `topic: versions` output. A warning is issued if no such topic is found.

## main\_nf\_version\_emit

The number of `topic: versions` outputs must
equal the number of `emit:` outputs whose name starts with `versions`.
A warning is issued if a legacy YAML-based `versions` emit is used instead
of a topic output.
