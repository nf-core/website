# files_exist

#### PipelineLint.files_exist()

Checks a given pipeline directory for required files.

Iterates through the pipeline’s directory content and checks that specified
files are either present or absent, as required.

#### NOTE

This test raises an `AssertionError` if neither `nextflow.config` or `main.nf` are found.
If these files are not found then this cannot be a Nextflow pipeline and something has gone badly wrong.
All lint tests are stopped immediately with a critical error message.

Files that _must_ be present:

```bash
.gitattributes
.gitignore
.nf-core.yml
.editorconfig
.prettierignore
.prettierrc.yml
.github/.dockstore.yml
.github/CONTRIBUTING.md
.github/ISSUE_TEMPLATE/bug_report.yml
.github/ISSUE_TEMPLATE/config.yml
.github/ISSUE_TEMPLATE/feature_request.yml
.github/PULL_REQUEST_TEMPLATE.md
.github/workflows/branch.yml
.github/workflows/ci.yml
.github/workflows/linting_comment.yml
.github/workflows/linting.yml
[LICENSE, LICENSE.md, LICENCE, LICENCE.md]  # NB: British / American spelling
assets/email_template.html
assets/email_template.txt
assets/nf-core-PIPELINE_logo_light.png
assets/sendmail_template.txt
conf/modules.config
conf/test.config
conf/test_full.config
CHANGELOG.md
CITATIONS.md
CODE_OF_CONDUCT.md
docs/images/nf-core-PIPELINE_logo_light.png
docs/images/nf-core-PIPELINE_logo_dark.png
docs/output.md
docs/README.md
docs/usage.md
lib/nfcore_external_java_deps.jar
lib/NfcoreTemplate.groovy
lib/Utils.groovy
lib/WorkflowMain.groovy
nextflow_schema.json
nextflow.config
README.md
```

Files that _should_ be present:

```bash
main.nf
assets/multiqc_config.yml
conf/base.config
conf/igenomes.config
.github/workflows/awstest.yml
.github/workflows/awsfulltest.yml
lib/WorkflowPIPELINE.groovy
pyproject.toml
```

Files that _must not_ be present, due to being renamed or removed in the template:

```bash
Singularity
parameters.settings.json
.nf-core.yaml  # NB: Should be yml, not yaml
bin/markdown_to_html.r
conf/aws.config
.github/workflows/push_dockerhub.yml
.github/ISSUE_TEMPLATE/bug_report.md
.github/ISSUE_TEMPLATE/feature_request.md
docs/images/nf-core-PIPELINE_logo.png
.markdownlint.yml
.yamllint.yml
lib/Checks.groovy
lib/Completion.groovy
lib/Workflow.groovy
```

Files that _should not_ be present:

```bash
.travis.yml
```
