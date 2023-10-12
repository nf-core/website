<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/files_unchanged.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.files_unchanged`

---

<a href="../../../../../../tools/nf_core/lint/files_unchanged.py#L12"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `files_unchanged`

```python
files_unchanged()
```

Checks that certain pipeline files are not modified from template output.

Iterates through the pipeline's directory content and compares specified files against output from the template using the pipeline's metadata. File content should not be modified / missing.

Files that must be unchanged:

```

     '.gitattributes',      '.github/.dockstore.yml',      '.github/CONTRIBUTING.md',      '.github/ISSUE_TEMPLATE/bug_report.md',      '.github/ISSUE_TEMPLATE/config.yml',      '.github/ISSUE_TEMPLATE/feature_request.md',      '.github/markdownlint.yml',      '.github/PULL_REQUEST_TEMPLATE.md',      '.github/workflows/branch.yml',      '.github/workflows/linting_comment.yml',      '.github/workflows/linting.yml',      'assets/email_template.html',      'assets/email_template.txt',      'assets/nf-core-PIPELINE_logo.png',      'assets/sendmail_template.txt',      'bin/markdown_to_html.py',      'CODE_OF_CONDUCT.md',      'docs/images/nf-core-PIPELINE_logo.png',      'docs/README.md',      'lib/nfcore_external_java_deps.jar'      'lib/NfcoreSchema.groovy',      ['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling

Files that can have additional content but must include the template contents:
```

     '.github/workflows/push_dockerhub_dev.yml',      '.github/workflows/push_dockerhub_release.yml',      '.gitignore',      'assets/multiqc_config.yaml',

```
.. tip:: You can configure the ``nf-core lint`` tests to ignore any of these checks by setting  the ``files_unchanged`` key as follows in your linting config file. For example:

 .. code-block:: yaml

 files_unchanged:
              - .github/workflows/branch.yml
              - assets/multiqc_config.yaml




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
```
