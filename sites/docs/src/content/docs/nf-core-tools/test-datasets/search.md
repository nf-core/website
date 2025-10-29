---
title: Search test datasets
subtitle: Search through module/pipeline test files in the nf-core/test-datasets repository on github
shortTitle: search
weight: 10
parentWeight: 30
---

The `nf-core test-datasets search` subcommand provides functionality to interactively search existing test-datasets for pipelines and modules from the commandline.
Test-datasets are hosted on the on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) github repository.

:::note
Not all files found on github can be searched with this subcommand.
Some auxiliary files like `README`, or `LICENSE` as well as all files starting with `.` are always ignored.
:::

## Search a branch for a data file

The subcommand fetches, filters and searches through the filetree of a branch of the test-dataset repository and prompts the user for a query, that can be tab-autocompleted.
A query can be specified as a command line argument to pre-populate the search field and if the query is unambiguous return the matching file without prompting.
A branch name to limit the search to is always required and if none is specified via the command line, the user will also be prompted to enter one.

![`nf-core test-datasets search --branch mag minigut_reads`](../../../../assets/images/tools/nf-core-test-datasets-search.svg)

:::note
To improve usability branch names can be entered via a tab-autocompletion. Alternatively, to list all branches see the [`list_branches` subcommand](/docs/nf-core-tools/test-datasets/list_branches).
:::

## Reusing the test-data files

As output from this subcommand, download links can be obtained or include statements for use in nf-core pipelines generated.
This uses the same flags `-u`/`--generate-dl-url` and `-p`/`--generate-nf-path` as the [`list` subcommand](/docs/nf-core-tools/test-datasets/list_branches).
