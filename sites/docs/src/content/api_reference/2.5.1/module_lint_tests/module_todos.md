# module_todos

#### `ModuleLint.module_todos(module){:python}`

Look for TODO statements in the module files

The nf-core module template contains a number of comment lines to help developers
of new modules know where they need to edit files and add content.
They typically have the following format:

```groovy
// TODO nf-core: Make some kind of change to the workflow here
```

..or in markdown:

```html
<!-- TODO nf-core: Add some detail to the docs here -->
```

This lint test runs through all files in the module and searches for these lines.
If any are found they will throw a warning.

:::note
Note that many GUI code editors have plugins to list all instances of _TODO_
in a given project directory. This is a very quick and convenient way to get
started on your pipeline!
:::
