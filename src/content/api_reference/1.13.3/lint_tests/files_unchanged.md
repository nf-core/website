# files_unchanged

#### `PipelineLint.files_unchanged(){:python}`

Checks that certain pipeline files are not modified from template output.

Iterates through the pipeline’s directory content and compares specified files
against output from the template using the pipeline’s metadata. File content
should not be modified / missing.

Files that must be unchanged:

```default
'.gitattributes',
'.github/.dockstore.yml',
'.github/CONTRIBUTING.md',
'.github/ISSUE_TEMPLATE/bug_report.md',
'.github/ISSUE_TEMPLATE/config.yml',
'.github/ISSUE_TEMPLATE/feature_request.md',
'.github/markdownlint.yml',
'.github/PULL_REQUEST_TEMPLATE.md',
'.github/workflows/branch.yml',
'.github/workflows/linting_comment.yml',
'.github/workflows/linting.yml',
'assets/email_template.html',
'assets/email_template.txt',
'assets/nf-core-PIPELINE_logo.png',
'assets/sendmail_template.txt',
'bin/markdown_to_html.py',
'CODE_OF_CONDUCT.md',
'docs/images/nf-core-PIPELINE_logo.png',
'docs/README.md',
'lib/nfcore_external_java_deps.jar'
'lib/NfcoreSchema.groovy',
['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling
```

Files that can have additional content but must include the template contents:

```default
'.github/workflows/push_dockerhub_dev.yml',
'.github/workflows/push_dockerhub_release.yml',
'.gitignore',
'assets/multiqc_config.yaml',
```

:::note
You can configure the `nf-core lint` tests to ignore any of these checks by setting
the `files_unchanged` key as follows in your linting config file. For example:

```yaml
files_unchanged:
  - .github/workflows/branch.yml
  - assets/multiqc_config.yaml
```

:::
