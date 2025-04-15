# readme

#### `PipelineLint.readme(){:python}`

Repository `README.md` tests

The `README.md` files for a project are very important and must meet some requirements:

- Nextflow badge
  - If no Nextflow badge is found, a warning is given
  - If a badge is found but the version doesn’t match the minimum version in the config file, the test fails
  - Example badge code:
    ```md
    [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
    ```

:::note
This badge are a markdown image `![alt-text](<image URL>)` _inside_ a markdown link `[markdown image](<link URL>)`, so a bit fiddly to write.
:::

- Zenodo release
  > - If pipeline is released but still contains a ‘zenodo.XXXXXXX’ tag, the test fails

To disable this test, add the following to the pipeline’s `.nf-core.yml` file:

To disable subsets of these tests, add the following to the pipeline’s `.nf-core.yml` file:

```yaml
lint:
  readme:
    - nextflow_badge
    - nfcore_template_badge
    - zenodo_release
```
