# actions_schema_validation

#### `PipelineLint.actions_schema_validation() → dict[str, list[str]]{:python}`

Checks that the GitHub Action workflow yml/yaml files adhere to the correct schema

nf-core pipelines use GitHub actions workflows to run CI tests, check formatting and also linting, among others.
These workflows are defined by `yml` scripts in `.github/workflows/`. This lint test verifies that these scripts are valid
by comparing them against the [JSON schema for GitHub workflows](https://www.schemastore.org/github-workflow).

To pass this test, make sure that all your workflows contain the required properties `on` and `jobs` and that
all other properties are of the correct type, as specified in the schema (link above).
