# nf_test_content

#### `PipelineLint.nf_test_content() → dict[str, list[str]]{:python}`

Checks that the pipeline nf-test files have the appropriate content.

This lint test checks the following files and content of these files:

- `*.nf.test` files should specify the `outdir` parameter:

```groovy
when {
    params {
        outdir = "$outputDir"
    }
}
```

- A versions.yml file should be included in the snapshot of all \*.nf.test files
- The nextflow.config file should contain:

  > ```groovy
  > modules_testdata_base_path = <path>
  > ```

  > and
  >
  > ```groovy
  > pipelines_testdata_base_path = <path>
  > ```

  > And should set the correct resource limits, as defined in the test profile

- The nf-test.config file should:
  : \* Make sure tests are relative to root directory <br/>
  ```groovy
  testsDir "."
  ```
    <br/>
    * Ensure a user-configurable nf-test directory
    <br/>
    ```groovy
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"
    ```
    <br/>
    * Use a test specific config
    <br/>
    ```groovy
    configFile "tests/nextflow.config"
    ```

All these checks can be skipped in the .nf-core.yml file using:

```yaml
lint:
  nf_test_content: False
```

or

```yaml
lint:
  nf_test_content:
    - tests/<test_name>.nf.test
    - tests/nextflow.config
    - nf-test.config
```
