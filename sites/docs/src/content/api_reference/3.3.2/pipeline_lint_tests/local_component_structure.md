# modules_structure

#### `PipelineLint.local_component_structure(){:python}`

Check that the local modules and subworkflows directories in a pipeline have the correct format:

```bash
modules/local/TOOL/SUBTOOL
```

Prior to nf-core/tools release 3.1.0 the directory structure allowed top-level \*.nf files:

```bash
modules/local/modules/TOOL_SUBTOOL.nf
```
