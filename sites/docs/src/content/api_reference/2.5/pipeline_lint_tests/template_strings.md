# template_strings

#### `PipelineLint.template_strings(){:python}`

Check for template placeholders.

The `nf-core create` pipeline template uses
[Jinja](https://jinja.palletsprojects.com/en/2.11.x/) behind the scenes.

This lint test fails if any Jinja template variables such as
`{{ pipeline_name }}` are found in your pipeline code.

Finding a placeholder like this means that something was probably copied and pasted
from the template without being properly rendered for your pipeline.

This test ignores any double-brackets prefixed with a dollar sign, such as
`${{ secrets.AWS_ACCESS_KEY_ID }}` as these placeholders are used in GitHub Actions workflows.
