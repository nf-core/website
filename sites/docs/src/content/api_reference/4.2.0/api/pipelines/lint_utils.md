# nf\_core.pipelines.utils

### `nf_core.pipelines.lint_utils.check_git_repo() → bool{:python}`

Check if the current directory is a git repository.

### `nf_core.pipelines.lint_utils.dump_json_with_prettier(file_name, file_content){:python}`

Dump a JSON file and run prettier on it.
:param file\_name: A file identifier as a string or pathlib.Path.
:type file\_name: Path | str
:param file\_content: Content to dump into the JSON file
:type file\_content: dict

### `nf_core.pipelines.lint_utils.dump_yaml_with_prettier(file_name: Path | str, file_content: dict) → None{:python}`

Dump a YAML file and run prettier on it.

- **Parameters:**
  - **file\_name** (_Path_ _|_ _str_) – A file identifier as a string or pathlib.Path.
  - **file\_content** (_dict_) – Content to dump into the YAML file

### `nf_core.pipelines.lint_utils.ignore_file(lint_name: str, file_path: Path, dir_path: Path) → list[list[str]]{:python}`

Ignore a file and add the result to the ignored list. Return the passed, failed, ignored and ignore\_configs lists.

### `nf_core.pipelines.lint_utils.print_fixes(lint_obj, plain_text=False){:python}`

Prints available and applied fixes

### `nf_core.pipelines.lint_utils.print_joint_summary(lint_obj, module_lint_obj, subworkflow_lint_obj, plain_text=False){:python}`

Print a joint summary of the general pipe lint tests and the module and subworkflow lint tests

### `nf_core.pipelines.lint_utils.print_results_plain_text(results_list, directory=None, component_type=None){:python}`

Print lint results in plain text format.

- **Parameters:**
  - **results\_list** – List of tuples (results, symbol, label, color, show\_condition)
  - **directory** – Base directory for relative paths (for component linting)
  - **component\_type** – “modules” or “subworkflows” (for component linting)

### `nf_core.pipelines.lint_utils.print_summary(rows, plain_text=False, summary_colour=None){:python}`

Print a summary table in plain text or rich format.

- **Parameters:**
  - **rows** – List of tuples (count, icon, label, color, always\_show)
  - **plain\_text** – If True, print in plain text format
  - **summary\_colour** – Color for the rich table border (default: auto based on failures)

### `nf_core.pipelines.lint_utils.run_prettier_on_file(file: Path | str | list[str]) → None{:python}`

Run the pre-commit hook prettier on a file.

- **Parameters:**
  **file** (_Path_ _|_ _str_) – A file identifier as a string or pathlib.Path.
- **Warns:**
  **If Prettier is not installed, a warning is logged.**
