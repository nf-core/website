---
title: 'nf-core/tools - 2.14.0'
subtitle: 'Spring cleaning :broom:'
pubDate: 2024-05-7T00:00:00+01:00
authors:
  - 'mirpedrol'
  - 'mashehu'
label:
  - 'tools'
---

This patch release contains some template changes and nf-core/tools updates. For a more detailed list of changes, you can read the [changelog](https://github.com/nf-core/tools/releases/tag/2.14.0).

# Highlights

## nf-core/tools functionalities

- We included a new linting test to assure that nf-test snapshots contain the versions.yml file.
- Components (modules and subworkflows) can now handle more possible git URLs (ssh:// and ftp://).

## Pipeline template

- We updated the GitHub Action which tests that the pipeline can be downloaded correctly (`download_pipeline.yml`):

  If you had errors with this GutHub test, they are fixed!

  First, the test will try to run a stub run of the downloaded pipeline. If this fails, because some modules don't have a stub test yet, it run the pipeline without the `-stub` option.

- We removed `pyproject.toml` from the template.

  This was used to lint Python code.
  The nf-core/pipeline template doesn't contain Python code anymore, thanks to the new utils subworkflows and the nf-validation (now nf-schema) plugin,
  which replace the Python script which was validating the input sample sheet.
  If you have other Python scripts in your pipeline and you would like to keep this linting, feel free to add this file back to your pipeline.

- Pipeline-specific institutional configs support is now activated for all pipelines by default.

- The `.nf-core.yml` file contains now the version of the pipeline template,
  corresponding to the version of nf-core/tools used for the last template update.

# How to merge the pipeline template updates

- Files inside the `.github` folder:

  These files are responsible for the Continuous Integration tests. In general, accept all changes made on these files.

  - `.github/workflows/ci.yml`: If you added your own tests to this file, for example, you added nf-test tests to your pipeline,
    keep your changes, but accept the template updates related to versions, e.g.

    ```diff "v1" "v2"
    - uses: nf-core/setup-nextflow@v1
    + uses: nf-core/setup-nextflow@v2
    ```

  - `.github/workflows/download_pipeline.yml`: This file is responsible for testing if the pipeline can be downloaded correctly. Accept the changes made to this file.

- `.nf-core-yml`:

  The version of nf-core used for the template update is added to the `.nf-core.yml` file.
  Accept the change of version.

- `README.md`

  We updated the link to Seqera Platform badge, accept this change.

- `assets/multiqc_config.yml`:

  Always accept changes made to this file before the line: `disable_version_detection: true`.
  Custom changes should be made after this line.

- `docs/usage.md`:

  We added a new profile `wave`. Accept this change.

- `nextflow.config`:

  We fixed a bug with `conda.channels`. Accept the changes made to this file.

- `pyproject.toml`:

  Python linting is now optional. If you have Python code on your pipeline and want to keep linting it, DON'T accept this change.
  Otherwise, it is safe to remove this file.

- Changes on `docs/`:

  Do NOT accept any change that removes custom docs that you added to your pipeline.

- Changes on CHANGELOG.md

  Do NOT accept any change which modified custom points of your `CHANGELOG.md`.

- `modules.json` and template modules and subworkflows:

  Do NOT accept any changes deleting your pipeline modules from `modules.json`.
  Template modules and subworkflows are updated on every template release. You can accept those changes and the changes on `modules.json` related to these.

  A safe way to add these changes is to NOT accept them.
  Then run `nf-core modules update` and `nf-core subworkflows update`.
  These commands will update all your modules and subworkflows and the `modules.json` file accordingly.
