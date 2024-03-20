---
title: Basic training to create an nf-core pipeline
subtitle: Working with Nextflow schema
---

## Nextflow Schema

All nf-core pipelines can be run with --help to see usage instructions. We can try this with the demo pipeline that we just created:

```

cd ../
nextflow run nf-core-demo/ --help

```

### Working with Nextflow schema

If you peek inside the nextflow_schema.json file you will see that it is quite an intimidating thing. The file is large and complex, and very easy to break if edited manually.

Thankfully, we provide a user-friendly tool for editing this file: nf-core schema build.

To see this in action, let’s add some new parameters to nextflow.config:

```

params {
demo = 'param-value-default'
foo = null
bar = false
baz = 12
// rest of the config file..

```

Then run nf-core schema build:

```

cd nf-core-demo/
nf-core schema build

```

The CLI tool should then prompt you to add each new parameter.
Here in the schema editor you can edit:

- Description and help text
- Type (string / boolean / integer etc)
- Grouping of parameters
- Whether a parameter is required, or hidden from help by default
- Enumerated values (choose from a list)
- Min / max values for numeric types
- Regular expressions for validation
- Special formats for strings, such as file-path
- Additional fields for files such as mime-type

:::tip{title="Exercise 7 - Using nextflow schema to add command line parameters"}

1.  Feed a string to your custom script from exercise 6 from the command line. Use `nf-core schema build` to add the parameter to the `nextflow.config` file.

      </details>

:::

<p class="text-center">
  <a href="/docs/contributing/nf_core_basic_training/add_custom_module/" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to Chapter 5
  </a>
  <a href="/docs/contributing/nf_core_basic_training/linting_modules/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 7 >
  </a>
</p>
