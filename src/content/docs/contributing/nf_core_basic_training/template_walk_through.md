---
title: Basic training to create an nf-core pipeline
subtitle: Exploring the nf-core template files
---

### Template code walk through

Now let us have a look at the files that were generated within the `nf-core-demotest` directory when we created this pipeline. You can see all files and directories either on the left hand side in the Explorer, or by running the command:

```bash
cd nf-core-demotest
tree
```

```
.
├── assets
│   ├── adaptivecard.json
│   ├── email_template.html
│   ├── email_template.txt
│   ├── methods_description_template.yml
│   ├── multiqc_config.yml
│   ├── nf-core-demotest_logo_light.png
│   ├── samplesheet.csv
│   ├── schema_input.json
│   ├── sendmail_template.txt
│   └── slackreport.json
├── bin
│   └── check_samplesheet.py
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── conf
│   ├── base.config
│   ├── igenomes.config
│   ├── modules.config
│   ├── test.config
│   └── test_full.config
├── docs
│   ├── images
│   │   ├── mqc_fastqc_adapter.png
│   │   ├── mqc_fastqc_counts.png
│   │   ├── mqc_fastqc_quality.png
│   │   ├── nf-core-demotest_logo_dark.png
│   │   └── nf-core-demotest_logo_light.png
│   ├── output.md
│   ├── README.md
│   └── usage.md
├── lib
│   ├── nfcore_external_java_deps.jar
│   ├── NfcoreTemplate.groovy
│   ├── Utils.groovy
│   ├── WorkflowDemotest.groovy
│   └── WorkflowMain.groovy
├── LICENSE
├── main.nf
├── modules
│   ├── local
│   │   └── samplesheet_check.nf
│   └── nf-core
│       ├── custom
│       │   └── dumpsoftwareversions
│       │       ├── main.nf
│       │       ├── meta.yml
│       │       └── templates
│       │           └── dumpsoftwareversions.py
│       ├── fastqc
│       │   ├── main.nf
│       │   └── meta.yml
│       └── multiqc
│           ├── main.nf
│           └── meta.yml
├── modules.json
├── nextflow.config
├── nextflow_schema.json
├── pyproject.toml
├── README.md
├── subworkflows
│   └── local
│       └── input_check.nf
├── tower.yml
└── workflows
    └── demotest.nf
```

These are the files in detail:

1. **main.nf**

   This file contains the main nextflow pipeline code. Mostly this file is not touched.

2. **workflows/demotest.nf**

   This file is where the pipeline is going to be assembled. It connects the different modules and subworkflows.

3. **CHANGELOG.md, CODE_OF_CONDUCT.md, LICENSE, README.md, CITATIONS.md**

   These are standard files created for github repositories. As a default, this pipeline will be under an MIT licence. The CODE_OF_CONDUCT is specific for the nf-core community.

4. **assets/**

   This directory contains different templates such as email templates or the MultiQC config. In this course contents of this directory can be largely ignored.

5. **bin/**

   The `bin` directory contains custom executable scripts, and is automatically added to the `PATH` by Nextflow allowing these scripts to become findable by Nextflow. As such, they can be called by name without using their absolute or relative path. The python script `check_samplesheet.py` is part of the nf-core template, since typically, nf-core pipelines require a samplesheet as one of their inputs.

6. **conf/** and **nextflow.config**

   The `nextflow.config` file is the main config file. In addition to supplying default parameters, it imports all the configurations in `conf/`. Importantly, `conf/` contains the `test.config` file, which is used for pipeline testing. In this course we are not going to touch config files, but they have been extensively covered in the following bytesize talks: [How nf-core configs work (nf-core/bytesize #2)](https://www.youtube.com/watch?v=cXBYusdjrc0&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=8&pp=gAQBiAQB), [Making a new institutional config profile (nf-core/bytesize #10)](https://www.youtube.com/watch?v=Ym1s6sKGzkw&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=9&pp=gAQBiAQB), [nf-core/bytesize: Using nf-core configs in custom pipelines](https://www.youtube.com/watch?v=zgcrI_0SUgg&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=40&pp=gAQBiAQB)

7. **docs/**

   This directory contains additional information to the README file. The most important files are the `output.md` and the `usage.md` files. `usage.md` should describe what exactly is needed to run the pipeline and `output.md` should describe all outputs that can be expected. Importantly, for nf-core pipelines, the information from these two files will automatically be displayed on the nf-core website page of the pipeline.

8. **lib/**

   This directory contains groovy functions and classes that are imported into the `main.nf` file to provide additional functionality not native to Nextflow.

9. **modules/local**

   This is where all your custom non-nf-core modules go. We will cover when and how to make local modules later in the course.

10. **modules/nf-core**

    All nf-core modules that are installed using the nf-core tooling will automatically show up in this directory. Keep them here, it is important for automatic updates.

11. **modules.json**

This file keeps track of modules installed using nf-core tools from the nf-core/modules repository. This file should only be updated using nf-core tools, and never manually.

12. **nextflow_schema.json**

    This file hosts all the parameters for the pipeline. Any new parameter should be added to this file using the `nf-core schema build` command. Similar to `output.md` and `usage.md`, the contents of `nextflow_schema.json` will get displayed on the pipeline page of nf-core pipelines.

13. **pyproject.toml**
14. **subworkflows/local**
15. **tower.yml**
16. **hidden directories and files**

    a) _.devcontainer/devcontainer.json_

    b) _.github/_

    Files in here are used for Continuous integration tests (CI) with github actions as well as other github related defaults, such as a template for issues. We will not touch on these in the course.

    c) _.gitignore_

    d) _.editorconfig_

    e) _.gitpod.yml_

    This file provides settings to create a Cloud development environment in your browser using Gitpod. It comes installed with the tools necessary to develop and test nf-core pipelines, modules, and subworkflows, allowing you to develop from anywhere without installing anything locally.

    f) _.nf-core.yml_

    g) _.pre-commit-config.yaml_

    h) _.prettierignore_

    i) _.prettierrc.yml_

:::tip{title="Exercise 2 - Test your knowledge of the nf-core pipeline structure"}

1. In which directory can you find the main script of the nf-core module `fastqc`
   <details>
      <summary>solution 1</summary>

   ```
   modules/nf-core/fastqc/
   ```

      </details>

2. Which file contains the main workflow of your new pipeline?
   <details>
      <summary>solution 2</summary>

   ```
   workflows/demotest.nf
   ```

      </details>

3. `check_samplesheet.py` is a script that can be called by any module of your pipeline, where is it located?
   <details>
      <summary>solution 3</summary>

   ```
   bin/
   ```

   This directory can also contain a custom scripts that you may wish to call from within a custom module.

      </details>

[MORE QUESTIONS CAN BE ADDED HERE]
:::

<p class="text-center">
  <a href="/docs/contributing/nf_core_basic_course/nf_core_create_tool/" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to Chapter 2
  </a>
  <a href="/docs/contributing/nf_core_basic_course/add_nf_core_module/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 4 >
  </a>
</p>
