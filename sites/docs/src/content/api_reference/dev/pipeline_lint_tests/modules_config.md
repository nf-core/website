# modules_config

#### `PipelineLint.modules_config() → Dict[str, List[str]]{:python}`

Make sure the conf/modules.config file follows the nf-core template, especially removed sections.

:::note
You can choose to ignore this lint tests by editing the file called
`.nf-core.yml` in the root of your pipeline and setting the test to false:

```yaml
lint:
  modules_config: False
```

:::

To disable this test only for specific modules, you can specify a list of module names.

```yaml
lint:
  modules_config:
    - fastqc
```
