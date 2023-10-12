<!-- markdownlint-disable -->

<a href="../../../../../../tools/nf_core/refgenie.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.refgenie`

Update a nextflow.config file with refgenie genomes

## **Global Variables**

- **NF_CFG_TEMPLATE**

---

<a href="../../../../../../tools/nf_core/refgenie.py#L104"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `update_config`

```python
update_config(rgc)
```

Update the genomes.config file after a local refgenie database has been updated

This function is executed after running 'refgenie pull <genome>/<asset>' The refgenie config file is transformed into a nextflow.config file, which is used to overwrited the 'refgenie_genomes.config' file. The path to the target config file is inferred from the following options, in order:

- the 'nextflow_config' attribute in the refgenie config file
- the NXF_REFGENIE_PATH environment variable
- otherwise defaults to: $NXF_HOME/nf-core/refgenie_genomes.config

Additionaly, an 'includeConfig' statement is added to the file $NXF_HOME/config

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
