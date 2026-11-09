---
title: Formatting
subtitle: Nextflow code formatting requirements
markdownPlugin: addNumbersToHeadings
weight: 9
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Nextflow syntax formatting

To maintain code quality and prevent issues, all code MUST be free of Nextflow warnings and errors.

Utilize the following command to check your code:

```bash
NXF_SYNTAX_PARSER=v2 nextflow lint modules/nf-core/module_name
```

:::note
The version two the of Nextflow strict syntax will be made default in Nextflow v26.04.
:::

Common issues within nf-core components to avoid:

| Old syntax                                                 | Prefered syntax                                                                |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------ |
| Unused `def args = task.ext.args ?: ''{:nextflow}` in stub | Delete the variable definition in the stub block                               |
| Undeclared variable `my_variable{:nextflow}`               | Specify `def` before the variable: `def my_variable{:nextflow}`                |
| `input.collect{ it[1].name }{:nextflow}`                   | Specify a variable name: `input.collect{ meta, file -> file.name }{:nextflow}` |
| `for{:nextflow}` loop                                      | Use the `.each` operator: `.each{}{:nextflow}` operator                        |

## General module code formatting

All code SHOULD otherwise be aligned to follow the automated Nextflow lint formatting (`nextflow lint -format -harshil-alignment`) or manually formatted following '[Harshil Alignment™️](/docs/contributing/code_editors_and_styling/harshil_alignment)', if it does not violate the Nextflow strict syntax specifications.
