# readme

#### PipelineLint.readme()

Repository `README.md` tests

The `README.md` files for a project are very important and must meet some requirements:

- Nextflow badge
  - If no Nextflow badge is found, a warning is given
  - If a badge is found but the version doesn’t match the minimum version in the config file, the test fails
  - Example badge code:
    ```md
    [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
    ```

#### NOTE

This badge are a markdown image `![alt-text](<image URL>)` _inside_ a markdown link `[markdown image](<link URL>)`, so a bit fiddly to write.

- Zenodo release
  > - If pipeline is released but still contains a ‘zenodo.XXXXXXX’ tag, the test fails
