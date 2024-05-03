# system_exit

#### `PipelineLint.system_exit(){:python}`

Check for System.exit calls in groovy/nextflow code

Calls to System.exit(1) should be replaced by throwing errors

This lint test looks for all calls to System.exit
in any file with the .nf or .groovy extension
