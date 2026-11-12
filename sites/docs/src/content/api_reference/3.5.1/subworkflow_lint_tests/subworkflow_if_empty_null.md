# subworkflow_if_empty_null

#### `SubworkflowLint.subworkflow_if_empty_null(subworkflow){:python}`

Check for ifEmpty(null)

There are two general cases for workflows to use the channel operator ifEmpty:
: 1. ifEmpty( \[ ] ) to ensure a process executes, for example when an input file is optional (although this can be replaced by toList()).
2\. When a channel should not be empty and throws an error ifEmpty { error … }, e.g. reading from an empty samplesheet.

There are multiple examples of workflows that inject null objects into channels using ifEmpty(null), which can cause unhandled null pointer exceptions.
This lint test throws warnings for those instances.
