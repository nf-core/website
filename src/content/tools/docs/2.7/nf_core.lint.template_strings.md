<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/lint/template_strings.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.lint.template_strings`

---

<a href="../../../../../../tools/nf_core/lint/template_strings.py#L6"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `template_strings`

```python
template_strings()
```

Check for template placeholders.

The `nf-core create` pipeline template uses `Jinja <https://jinja.palletsprojects.com/en/2.11.x/>`\_ behind the scenes.

This lint test fails if any Jinja template variables such as `{{ pipeline_name }}` are found in your pipeline code.

Finding a placeholder like this means that something was probably copied and pasted from the template without being properly rendered for your pipeline.

This test ignores any double-brackets prefixed with a dollar sign, such as `${{ secrets.AWS_ACCESS_KEY_ID }}` as these placeholders are used in GitHub Actions workflows.

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
