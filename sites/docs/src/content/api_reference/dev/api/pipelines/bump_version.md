# nf_core.pipelines.bump_version

Bumps the version number in all appropriate files for
a nf-core pipeline.

### `nf_core.pipelines.bump_version.bump_nextflow_version(pipeline_obj:{:python}`[`Pipeline{:python}`](../utils#nf_core.utils.Pipeline)`, new_version: str) → None{:python}`

Bumps the required Nextflow version number of a pipeline.

- **Parameters:**
  - **pipeline_obj** ([_nf_core.utils.Pipeline_](../utils#nf_core.utils.Pipeline)) – A Pipeline object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the required Nextflow version.

### `nf_core.pipelines.bump_version.bump_pipeline_version(pipeline_obj:{:python}`[`Pipeline{:python}`](../utils#nf_core.utils.Pipeline)`, new_version: str) → None{:python}`

Bumps a pipeline version number.

- **Parameters:**
  - **pipeline_obj** ([_nf_core.utils.Pipeline_](../utils#nf_core.utils.Pipeline)) – A Pipeline object that holds information
    about the pipeline contents and build files.
  - **new_version** (_str_) – The new version tag for the pipeline. Semantic versioning only.

### `nf_core.pipelines.bump_version.handle_error(message: str, required: bool){:python}`

### `nf_core.pipelines.bump_version.log_change(old_content: str, new_content: str){:python}`

### `nf_core.pipelines.bump_version.update_file_version(filename: str | Path, pipeline_obj:{:python}`[`Pipeline{:python}`](../utils#nf_core.utils.Pipeline)`, patterns: list[tuple[str, str]], required: bool = True, yaml_key: list[str] | None = None) → None{:python}`

Updates a file with a new version number.

- **Parameters:**
  - **filename** (_str_) – The name of the file to update.
  - **pipeline_obj** ([_nf_core.utils.Pipeline_](../utils#nf_core.utils.Pipeline)) – A Pipeline object that holds information
    about the pipeline contents.
  - **patterns** (_List_ \*\[\*_Tuple_ \*\[\*_str_ _,_ _str_ _]_ _]_) – A list of tuples containing the regex patterns to
    match and the replacement strings.
  - **required** (_bool_ _,_ _optional_) – Whether the file is required to exist. Defaults to True.
  - **yaml_key** (_List_ \*\[\*_str_ _]_ _|_ _None_ _,_ _optional_) – The YAML key to update. Defaults to None.

### `nf_core.pipelines.bump_version.update_text_file(fn: Path, patterns: list[tuple[str, str]], required: bool){:python}`

Updates a text file with a new version number.

- **Parameters:**
  - **fn** (_Path_) – The name of the file to update.
  - **patterns** (_List_ \*\[\*_Tuple_ \*\[\*_str_ _,_ _str_ _]_ _]_) – A list of tuples containing the regex patterns to
    match and the replacement strings.
  - **required** (_bool_) – Whether the file is required to exist.

### `nf_core.pipelines.bump_version.update_yaml_file(fn: Path, patterns: list[tuple[str, str]], yaml_key: list[str], required: bool){:python}`

Updates a YAML file with a new version number.

- **Parameters:**
  - **fn** (_Path_) – The name of the file to update.
  - **patterns** (_List_ \*\[\*_Tuple_ \*\[\*_str_ _,_ _str_ _]_ _]_) – A list of tuples containing the regex patterns to
    match and the replacement strings.
  - **yaml_key** (_List_ \*\[\*_str_ _]_) – The YAML key to update.
  - **required** (_bool_) – Whether the file is required to exist.
