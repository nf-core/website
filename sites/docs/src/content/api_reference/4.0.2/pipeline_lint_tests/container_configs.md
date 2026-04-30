# container_configs

#### `PipelineLint.container_configs(){:python}`

Check that the container configuration files in `conf/` are up to date.

Runs `nextflow inspect` to regenerate container configuration files directly
in `conf/` and uses `git diff` to detect changes. If not in `--fix` mode
the working tree is restored to its original state afterwards.

Can be skipped by adding the following to the `.nf-core.yml` file:

```yaml
lint:
  container_configs: False
```
