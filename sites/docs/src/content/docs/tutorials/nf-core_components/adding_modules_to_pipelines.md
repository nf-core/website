---
title: "Adding Modules to Pipelines"
subtitle: Quick reference of the steps needed to add a new module to a pipeline.
weight: 20
---

This tutorial is a quick reference overview of the steps you will normally take when adding a new module to a pipeline.

Different pipelines may have different workflows, but the following steps will cover the most common aspects.

1. Install module to the pipeline

   ```bash
   nf-core modules install <tool>/<subtool>
   ```

2. Add `include` statement at the top of the pipeline's (sub)workflow to import module (as suggested by `nf-core modules install`)
3. Insert module execution into the relevant place (sub)workflow script `TOOL(ch_input)`
   - Make sure to mix in the module's version into the version channel, e.g. `ch_versions = ch_versions.mix(TOOL.out.versions)`
   - Make sure to mix any output files for MultiQC into a relevant channel, e.g. `ch_multiqc_files = ch_multiqc.mix(TOOL.out.log)`
4. Create a section in `conf/modules.conf` for the module (with a default results directory and output file pattern)
5. Add any necessary parameters for the module with defaults to `nextflow.config`
6. Insert any pipeline level parameters (`params.*`) into the `ext.args` of corresponding `conf/modules.config`
   - In some cases these may need to be passed directly to the module itself, e.g. `FASTP( reads, params.save_trimmed_fail, params.save_merged)`
   - Ensure all directives which use `params.*` supply the value as a closure (i.e., enclosed within `{}`, e.g., `ext.args = { "--option $params.option" }`).
     This ensures parameters supplied in a config `-c` are correctly resolved.
7. Update the `nextflow_schema.json` to include new parameters with

   ```bash
   nf-core schema build
   ```

8. Update `assets/multiqc_config.yml` to include any new MultiQC modules (if any exist) and specify order in the report
9. Add a citation for the new tool/module to `citations.md`
10. Update `docs/USAGE.md` to describe any important information about running of the module (this can be optional in some cases)
11. Update `docs/OUTPUT.md` to describe the directories output files of the module
12. Update `README.md` mentioning the tool is used and any pipeline diagrams (optional)
13. If not already installed, [install](https://prettier.io) [prettier](https://nf-co.re/events/2022/bytesize-41-prettier) (prettier can also be installed using Conda) and then run it formatting on the whole repository

    ```bash
    prettier -w .
    ```

    - Or you can run [pre-commit](https://nf-co.re/events/2023/bytesize_precommit) to use prettier when committing your code

    - If you forget this step you can also post a comment on the open PR once made with `@nf-core-bot fix linting`

14. Run a local test of the pipeline with the included new functionality to check it works.

    ```bash
    mkdir test/ && cd test/
    nextflow run ../main.nf -profile test,<docker,singularity,conda> --outdir ./results <include new parameters required to activate new functionality if necessary>
    ```

15. Lint the new code with

    ```bash
    nf-core lint
    ```

Then open the pull request!
