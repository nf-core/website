---
title: 'Tutorial: Adding Modules to Pipelines'
subtitle: Quick reference of the steps needed to add a new module to a pipeline.
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
7. Update the `nextflow_schema.json` to include new parameters with

   ```bash
   nf-core schema build
   ```

8. Update `assets/multiqc_config.yml` to include any new MultiQC modules (if any exist) and specify order in the report
9. Add a citation for the new tool/module to `citations.md`
10. Update `docs/USAGE.md` to describe any important information about running of the module (this can be optional in some cases)
11. Update `docs/OUTPUT.md` to describe the directories output files of the module
12. Update `README.md` mentioning the tool is used and any pipeline diagrams (optional)
13. Run prettier formatting on the whole repository (if installed)

    ```bash
    prettier -w .
    ```

    - Or run VSCode 'Format Document' function on: `citations.md`, `usage.md`, `output.md`, `README.md`, `modules.json`, `nextflow_schema.json`
    - If you forget this step you can also post a comment on the open PR once made with `@nf-core-bot fix linting`

14. Lint the new code with

    ```bash
    nf-core lint
    ```

Then open the pull request!
