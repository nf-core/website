# base_config

#### `PipelineLint.base_config() → dict[str, list[str]]{:python}`

Make sure the conf/base.config file follows the nf-core template, especially removed sections.

:::note
You can choose to ignore this lint tests by editing the file called
`.nf-core.yml` in the root of your pipeline and setting the test to false:

```yaml
lint:
  base_config: False
```

:::
