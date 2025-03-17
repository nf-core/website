# merge_markers

#### `PipelineLint.merge_markers(){:python}`

Check for remaining merge markers.

This test looks for remaining merge markers in the code, e.g.:
`>>>>>>>` or `<<<<<<<`

:::note
You can choose to ignore this lint tests by editing the file called
`.nf-core.yml` in the root of your pipeline and setting the test to false:

```yaml
lint:
  merge_markers: False
```

:::

To disable this test only for specific files, you can specify a list of file paths to ignore.
For example, to ignore a pdf you added to the docs:

```yaml
lint:
  merge_markers:
    - docs/my_pdf.pdf
```
