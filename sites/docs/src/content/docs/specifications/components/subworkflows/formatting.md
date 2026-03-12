---
title: Formatting
subtitle: Nextflow code formatting requirements
markdownPlugin: addNumbersToHeadings
weight: 7
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

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

| Old syntax                                 | Prefered syntax                                                              |
| ------------------------------------------ | ---------------------------------------------------------------------------- |
| Undeclared variable `my_variable{:groovy}` | Specify `def` before the variable: `def my_variable{:groovy}`                |
| `input.collect{ it[1].name }{:groovy}`     | Specify a variable name: `input.collect{ meta, file -> file.name }{:groovy}` |
| `for{:groovy}` loop                        | Use the `.each` operator: `.each{}{:groovy}` operator                        |

## General module code formatting

All code SHOULD otherwise be aligned to follow the '[Harshil Alignment™️](/docs/contributing/code_editors_and_styling/harshil_alignment)' format, if it does not violate the Nextflow strict syntax specifications.
