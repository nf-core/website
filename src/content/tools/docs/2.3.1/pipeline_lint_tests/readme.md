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
- Bioconda badge
  - If your pipeline contains a file called `environment.yml` in the root directory, a bioconda badge is required
  - Required badge code:
    ```md
    [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
    ```

#### NOTE

These badges are a markdown image `![alt-text](<image URL>)` _inside_ a markdown link `[markdown image](<link URL>)`, so a bit fiddly to write.
