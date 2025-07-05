# version_consistency

#### `PipelineLint.version_consistency(){:python}`

Pipeline and container version number consistency.

:::note
This test only runs when the `--release` flag is set for `nf-core pipelines lint`,
or `$GITHUB_REF` is equal to `main`.
:::

This lint fetches the pipeline version number from four possible locations:

- The pipeline config, `manifest.version`
- The docker container in the pipeline config, `process.container`
  > Some pipelines may not have this set on a pipeline level. If it is not found, it is ignored.
- `$GITHUB_REF`, if it looks like a release tag (`refs/tags/<something>`)
- The YAML file .nf-core.yml

The test then checks that:

- The container name has a tag specified (eg. `nfcore/pipeline:version`)
- The pipeline version number is numeric (contains only numbers and dots)
- That the version numbers all match one another
