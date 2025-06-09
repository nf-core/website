# files_unchanged

#### `PipelineLint.files_unchanged() → Dict[str, List[str] | bool]{:python}`

Checks that certain pipeline files are not modified from template output.

Iterates through the pipeline’s directory content and compares specified files
against output from the template using the pipeline’s metadata. File content
should not be modified / missing.

Files that must be unchanged:

```default
.gitattributes
.prettierrc.yml
.github/.dockstore.yml
.github/CONTRIBUTING.md
.github/ISSUE_TEMPLATE/bug_report.yml
.github/ISSUE_TEMPLATE/config.yml
.github/ISSUE_TEMPLATE/feature_request.yml
.github/PULL_REQUEST_TEMPLATE.md
.github/workflows/branch.yml
.github/workflows/linting_comment.yml
.github/workflows/linting.yml
assets/email_template.html
assets/email_template.txt
assets/nf-core-PIPELINE_logo_light.png
assets/sendmail_template.txt
CODE_OF_CONDUCT.md
docs/images/nf-core-PIPELINE_logo_light.png
docs/images/nf-core-PIPELINE_logo_dark.png
docs/README.md'
['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling
```

Files that can have additional content but must include the template contents:

```default
.gitignore
.prettierignore
```

:::note
You can configure the `nf-core pipelines lint` tests to ignore any of these checks by setting
the `files_unchanged` key as follows in your `.nf-core.yml` config file. For example:

```yaml
lint:
  files_unchanged:
    - .github/workflows/branch.yml
```

:::
