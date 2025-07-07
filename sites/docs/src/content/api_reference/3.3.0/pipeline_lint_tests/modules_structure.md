# modules_structure

#### `PipelineLint.modules_structure(){:python}`

Check that the structure of the modules directory in a pipeline is the correct one:

```bash
modules/nf-core/TOOL/SUBTOOL
```

Prior to nf-core/tools release 2.6 the directory structure had an additional level of nesting:

```bash
modules/nf-core/modules/TOOL/SUBTOOL
```
