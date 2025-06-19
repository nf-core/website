# files_exist

#### `PipelineLint.files_exist() → dict[str, list[str]]{:python}`

Checks a given pipeline directory for required files.

Iterates through the pipeline’s directory content and checks that specified
files are either present or absent, as required.

:::note
This test raises an `AssertionError` if neither `nextflow.config` or `main.nf` are found.
If these files are not found then this cannot be a Nextflow pipeline and something has gone badly wrong.
All lint tests are stopped immediately with a critical error message.
:::

Files that _must_ be present:

```bash
.gitattributes
.gitignore
.nf-core.yml
.prettierignore
.prettierrc.yml
.github/.dockstore.yml
.github/CONTRIBUTING.md
.github/ISSUE_TEMPLATE/bug_report.yml
.github/ISSUE_TEMPLATE/config.yml
.github/ISSUE_TEMPLATE/feature_request.yml
.github/PULL_REQUEST_TEMPLATE.md
.github/workflows/branch.yml
.github/workflows/nf-test.yml
.github/actions/get-shards/action.yml
.github/actions/nf-test/action.yml
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
nextflow_schema.json
nextflow.config
nf-test.config
README.md
tests/default.nf.test
```

Files that _should_ be present:

```bash
main.nf
assets/multiqc_config.yml
conf/base.config
conf/igenomes.config
.github/workflows/awstest.yml
.github/workflows/awsfulltest.yml
ro-crate-metadata.json
```

Files that _must not_ be present, due to being renamed or removed in the template:

```bash
.github/ISSUE_TEMPLATE/bug_report.md
.github/ISSUE_TEMPLATE/feature_request.md
.github/workflows/push_dockerhub.yml
.markdownlint.yml
.nf-core.yaml  # NB: Should be yml, not yaml
.yamllint.yml
bin/markdown_to_html.r
conf/aws.config
docs/images/nf-core-PIPELINE_logo.png
lib/Checks.groovy
lib/Completion.groovy
lib/NfcoreTemplate.groovy
lib/Utils.groovy
lib/Workflow.groovy
lib/WorkflowMain.groovy
lib/WorkflowPIPELINE.groovy
lib/nfcore_external_java_deps.jar
parameters.settings.json
pipeline_template.yml # saving information in .nf-core.yml
Singularity
```

Files that _should not_ be present:

```bash
.travis.yml
```

:::note
You can configure the `nf-core pipelines lint` tests to ignore any of these checks by setting
the `files_exist` key as follows in your `.nf-core.yml` config file. For example:

```yaml

```

:::

lint:
: files_exist:
: - assets/multiqc_config.yml
