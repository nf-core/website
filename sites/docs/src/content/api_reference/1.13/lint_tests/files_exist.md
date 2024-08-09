# files_exist

#### `PipelineLint.files_exist(){:python}`

Checks a given pipeline directory for required files.

Iterates through the pipeline’s directory content and checks that specified
files are either present or absent, as required.

:::note
This test raises an `AssertionError` if neither `nextflow.config` or `main.nf` are found.
If these files are not found then this cannot be a Nextflow pipeline and something has gone badly wrong.
All lint tests are stopped immediately with a critical error message.
:::

Files that **must** be present:

```default
.gitattributes
.github/.dockstore.yml
.github/CONTRIBUTING.md
.github/ISSUE_TEMPLATE/bug_report.md
.github/ISSUE_TEMPLATE/config.yml
.github/ISSUE_TEMPLATE/feature_request.md
.github/markdownlint.yml
.github/PULL_REQUEST_TEMPLATE.md
.github/workflows/branch.yml
.github/workflows/ci.yml
.github/workflows/linting_comment.yml
.github/workflows/linting.yml
[LICENSE, LICENSE.md, LICENCE, LICENCE.md]  # NB: British / American spelling
assets/email_template.html
assets/email_template.txt
assets/nf-core-PIPELINE_logo.png
assets/sendmail_template.txt
bin/markdown_to_html.py
CHANGELOG.md
CODE_OF_CONDUCT.md
CODE_OF_CONDUCT.md
docs/images/nf-core-PIPELINE_logo.png
docs/output.md
docs/README.md
docs/README.md
docs/usage.md
lib/nfcore_external_java_deps.jar
lib/NfcoreSchema.groovy
nextflow_schema.json
nextflow.config
README.md
```

Files that _should_ be present:

```default
'main.nf',
'environment.yml',
'Dockerfile',
'conf/base.config',
'.github/workflows/awstest.yml',
'.github/workflows/awsfulltest.yml'
```

Files that _must not_ be present:

```default
'Singularity',
'parameters.settings.json',
'bin/markdown_to_html.r',
'conf/aws.config',
'.github/workflows/push_dockerhub.yml'
```

Files that _should not_ be present:

```default
'.travis.yml'
```
