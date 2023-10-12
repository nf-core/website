<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/files_unchanged.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.files_unchanged`

---

<a href="../../../../../../tools/nf_core/lint/files_unchanged.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `files_unchanged`

```python
files_unchanged()
```

Checks that certain pipeline files are not modified from template output.

Iterates through the pipeline's directory content and compares specified files against output from the template using the pipeline's metadata. File content should not be modified / missing.

Files that must be unchanged:

```

     .gitattributes      .prettierrc.yml      .github/.dockstore.yml      .github/CONTRIBUTING.md      .github/ISSUE_TEMPLATE/bug_report.yml      .github/ISSUE_TEMPLATE/config.yml      .github/ISSUE_TEMPLATE/feature_request.yml      .github/PULL_REQUEST_TEMPLATE.md      .github/workflows/branch.yml      .github/workflows/linting_comment.yml      .github/workflows/linting.yml      assets/email_template.html      assets/email_template.txt      assets/nf-core-PIPELINE_logo_light.png      assets/sendmail_template.txt      CODE_OF_CONDUCT.md      docs/images/nf-core-PIPELINE_logo_light.png      docs/images/nf-core-PIPELINE_logo_dark.png      docs/README.md'      lib/nfcore_external_java_deps.jar      lib/NfcoreSchema.groovy      lib/NfcoreTemplate.groovy      ['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling

Files that can have additional content but must include the template contents:
```

     .gitignore      .prettierignore

```
.. tip:: You can configure the ``nf-core lint`` tests to ignore any of these checks by setting  the ``files_unchanged`` key as follows in your ``.nf-core.yml`` config file. For example:

 .. code-block:: yaml

 lint:  files_unchanged:
                    - .github/workflows/branch.yml




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
```
