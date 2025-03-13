---
title: Test a module
subtitle: Create a test for a module
shortTitle: test
weight: 90
---

All modules on [nf-core/modules](https://github.com/nf-core/modules) have a strict requirement of being unit tested using minimal test data. We use [nf-test](https://code.askimed.com/nf-test/) as our testing framework.
Each module comes already with a template for the test file in `test/main.nf.test`. Replace the placeholder code in that file with your specific input, output and proces. In order to generate the corresponding snapshot after writing your test, you can use the `nf-core modules test` command. This command will run `nf-test test` twice, to also check for snapshot stability, i.e. that the same snapshot is generated on multiple runs.

You can specify the module name in the form TOOL/SUBTOOL in the command or provide it later through interactive prompts.

<!-- RICH-CODEX
working_dir: tmp/modules
timeout: 30
extra_env:
  PROFILE: 'conda'
-->

![`nf-core modules test fastqc --no-prompts`](/images/tools/nf-core-modules-test.svg)

In case you changed something in the test and want to update the snapshot, run

```bash
nf-core modules test --update
```

If you want to run the test only once without checking for snapshot stability, you can use the `--once` flag.
