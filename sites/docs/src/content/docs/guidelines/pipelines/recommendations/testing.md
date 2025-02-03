---
title: Testing
subtitle: Use nf-test to validate pipeline
weight: 210
---

Pipelines should use nf-test for valid pipeline testing.

All pipelines should support nf-test and use it for pipeline level tests. Additional tests for local subworkflows and modules are recommended. Modules and subworkflows from nf-core/modules should include nf-tests which can also be used within the pipeline.

Pipeline level tests can facilitate more reliable and reproducible pipelines by ensuring the pipeline produces identical results with every run. You must add pipeline tests that work with `-profile test` and you should reuse this profile within one nf-test.

### Pipeline nf-test overview and structure

Within the base directory of the repository, there is a configuration file for nf-test, named `nf-test.config`. This will set the following options:

- Set the `testsDir` to the base of the repository so it includes all files
- Set the default profile(s) for nf-test to include `test` (this can be overridden on the command line)
- Add an additional configuration file specific for nf-test located in `tests/nextflow.config`

Within the nf-test specific configuration file, you can add specific requirements for running nf-tests but this should not include parameters or options as these should be available in all contexts.

Pipeline level tests are located in the `tests` directory. Each nf-test file must contain a single pipeline test which tests for a single scenario. Although nf-test supports multiple tests within a single nf-test file, keeping thme in separate files makes it easier to launch individual pipeline tests. Each nf-test file should be named after the scenario it tests in the following format:

```tree
tests
├── single_end.main.nf.test
└── paired_end.main.nf.test
```

### Pipeline nf-tests additional guidance

The same guidelines for test profiles, test data and nf-test also apply to pipeline tests. In addition, the following guidelines apply:

- To ensure all output files are caught, the `params.outdir` should be set the the nf-test variable `outputDir`
- The tag `PIPELINE` and the pipeline name should be added to all tests
