# container\_configs

#### `PipelineLint.container_configs(){:python}`

Check that the container configuration files in `conf/` are up to date.

Scans all `meta.yml` files under `modules/` that contain a `containers`
key, reads the process name from the sibling `main.nf`, and regenerates
the container configuration files in `conf/`. Uses `git diff` to detect
changes. If not in `--fix` mode the working tree is restored to its
original state afterwards.

Can be skipped by adding the following to the `.nf-core.yml` file:

```yaml
lint:
  container_configs: False
```
