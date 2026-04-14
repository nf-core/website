# main_nf

#### `ModuleLint.main_nf(module: NFCoreComponent, fix_version: bool, registry: str, progress_bar: Progress) → tuple[list[str], list[str]]{:python}`

Lint a `main.nf` module file

Can also be used to lint local module files,
in which case failures will be reported as warnings.

The following checks are performed:

- `main_nf_exists`: The `main.nf` file must exist.
- `deprecated_dsl2`: The file must not contain deprecated DSL2 identifiers
  (`initOptions`, `saveFiles`, `getSoftwareName`, `getProcessName`,
  `publishDir`).
- `main_nf_script_outputs`: The process must have an `output:` block.
- `main_nf_container`: Container tags across the `singularity`, `docker`,
  and `conda` directives must reference the same software version. A warning
  is issued if they do not match.
- `main_nf_script_shell`: Exactly one of `script:`, `shell:`, or `exec:`
  blocks must be present.
- `main_nf_shell_template`: If a `shell:` block is used, it must call
  a `template`.
- `main_nf_meta_output`: If `meta` is present in the module inputs, it
  must also appear in at least one output channel.
- `main_nf_version_topic`: The module should emit software versions using
  a `topic: versions` output. A warning is issued if no such topic is found.
- `main_nf_version_emit`: The number of `topic: versions` outputs must
  equal the number of `emit:` outputs whose name starts with `versions`.
  A warning is issued if a legacy YAML-based `versions` emit is used instead
  of a topic output.
