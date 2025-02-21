# nf_core.utils

Common utility functions for the nf-core python package.

### _`class{:python}`_`nf_core.utils.GitHubAPISession{:python}`

Bases: `CachedSession`

Class to provide a single session for interacting with the GitHub API for a run.
Inherits the requests_cache.CachedSession and adds additional functionality,
such as automatically setting up GitHub authentication if we can.

#### `get(url, **kwargs){:python}`

Initialise the session if we haven’t already, then call the superclass get method.

#### `lazy_init() → None{:python}`

Initialise the object.

Only do this when it’s actually being used (due to global import)

#### `log_content_headers(request, post_data=None){:python}`

Try to dump everything to the console, useful when things go wrong.

#### `request_retry(url, post_data=None){:python}`

Try to fetch a URL, keep retrying if we get a certain return code.

Used in nf-core pipelines sync code because we get 403 errors: too many simultaneous requests
See <https://github.com/nf-core/tools/issues/911>

#### `safe_get(url){:python}`

Run a GET request, raise a nice exception with lots of logging if it fails.

#### `setup_github_auth(auth=None){:python}`

Try to automatically set up GitHub authentication

### _`pydantic model{:python}`_`nf_core.utils.NFCoreTemplateConfig{:python}`

Bases: `BaseModel`

Template configuration schema

<p><details  class="autodoc_pydantic_collapsable_json">
<summary>Show JSON schema</summary>
```json
{
   "title": "NFCoreTemplateConfig",
   "description": "Template configuration schema",
   "type": "object",
   "properties": {
      "org": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Org"
      },
      "name": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Name"
      },
      "description": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Description"
      },
      "author": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Author"
      },
      "version": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Version"
      },
      "force": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": true,
         "title": "Force"
      },
      "outdir": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "format": "path",
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Outdir"
      },
      "skip_features": {
         "anyOf": [
            {
               "items": {},
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Skip Features"
      },
      "is_nfcore": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Is Nfcore"
      }
   }
}
```

</details></p>
* **Fields:**
  - [`author (str | None)`](#nf_core.utils.NFCoreTemplateConfig.author)
  - [`description (str | None)`](#nf_core.utils.NFCoreTemplateConfig.description)
  - [`force (bool | None)`](#nf_core.utils.NFCoreTemplateConfig.force)
  - [`is_nfcore (bool | None)`](#nf_core.utils.NFCoreTemplateConfig.is_nfcore)
  - [`name (str | None)`](#nf_core.utils.NFCoreTemplateConfig.name)
  - [`org (str | None)`](#nf_core.utils.NFCoreTemplateConfig.org)
  - [`outdir (str | pathlib.Path | None)`](#nf_core.utils.NFCoreTemplateConfig.outdir)
  - [`skip_features (list | None)`](#nf_core.utils.NFCoreTemplateConfig.skip_features)
  - [`version (str | None)`](#nf_core.utils.NFCoreTemplateConfig.version)
* **Validators:**
  - [`outdir_to_str`](#nf_core.utils.NFCoreTemplateConfig.outdir_to_str) » [`outdir`](#nf_core.utils.NFCoreTemplateConfig.outdir)

#### _`field{:python}`_`author{:python}`_: str | None_`{:python}`_= None_

Pipeline author

#### _`field{:python}`_`description{:python}`_: str | None_`{:python}`_= None_

Pipeline description

#### _`field{:python}`_`force{:python}`_: bool | None_`{:python}`_= True_

Force overwrite of existing files

#### _`field{:python}`_`is_nfcore{:python}`_: bool | None_`{:python}`_= None_

Whether the pipeline is an nf-core pipeline.

#### _`field{:python}`_`name{:python}`_: str | None_`{:python}`_= None_

Pipeline name

#### _`field{:python}`_`org{:python}`_: str | None_`{:python}`_= None_

Organisation name

#### _`field{:python}`_`outdir{:python}`_: str | Path | None_`{:python}`_= None_

Output directory

- **Validated by:**
  - [`outdir_to_str`](#nf_core.utils.NFCoreTemplateConfig.outdir_to_str)

#### _`field{:python}`_`skip_features{:python}`_: list | None_`{:python}`_= None_

Skip features. See <https://nf-co.re/docs/nf-core-tools/pipelines/create> for a list of features.

#### _`field{:python}`_`version{:python}`_: str | None_`{:python}`_= None_

Pipeline version

#### _`validator{:python}`_`outdir_to_str{:python}`_»_`{:python}`[_outdir_](#nf_core.utils.NFCoreTemplateConfig.outdir)

#### `get(item: str, default: Any = None) → Any{:python}`

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

### _`pydantic model{:python}`_`nf_core.utils.NFCoreYamlConfig{:python}`

Bases: `BaseModel`

.nf-core.yml configuration file schema

<p><details  class="autodoc_pydantic_collapsable_json">
<summary>Show JSON schema</summary>
```json
{
   "title": "NFCoreYamlConfig",
   "description": ".nf-core.yml configuration file schema",
   "type": "object",
   "properties": {
      "repository_type": {
         "anyOf": [
            {
               "enum": [
                  "pipeline",
                  "modules"
               ],
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Repository Type"
      },
      "nf_core_version": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Nf Core Version"
      },
      "org_path": {
         "anyOf": [
            {
               "type": "string"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Org Path"
      },
      "lint": {
         "anyOf": [
            {
               "$ref": "#/$defs/NFCoreYamlLintConfig"
            },
            {
               "type": "null"
            }
         ],
         "default": null
      },
      "template": {
         "anyOf": [
            {
               "$ref": "#/$defs/NFCoreTemplateConfig"
            },
            {
               "type": "null"
            }
         ],
         "default": null
      },
      "bump_version": {
         "anyOf": [
            {
               "additionalProperties": {
                  "type": "boolean"
               },
               "type": "object"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Bump Version"
      },
      "update": {
         "anyOf": [
            {
               "additionalProperties": {
                  "anyOf": [
                     {
                        "type": "string"
                     },
                     {
                        "type": "boolean"
                     },
                     {
                        "additionalProperties": {
                           "anyOf": [
                              {
                                 "type": "string"
                              },
                              {
                                 "additionalProperties": {
                                    "anyOf": [
                                       {
                                          "type": "string"
                                       },
                                       {
                                          "type": "boolean"
                                       }
                                    ]
                                 },
                                 "type": "object"
                              }
                           ]
                        },
                        "type": "object"
                     }
                  ]
               },
               "type": "object"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Update"
      }
   },
   "$defs": {
      "NFCoreTemplateConfig": {
         "description": "Template configuration schema",
         "properties": {
            "org": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Org"
            },
            "name": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Name"
            },
            "description": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Description"
            },
            "author": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Author"
            },
            "version": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Version"
            },
            "force": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": true,
               "title": "Force"
            },
            "outdir": {
               "anyOf": [
                  {
                     "type": "string"
                  },
                  {
                     "format": "path",
                     "type": "string"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Outdir"
            },
            "skip_features": {
               "anyOf": [
                  {
                     "items": {},
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Skip Features"
            },
            "is_nfcore": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Is Nfcore"
            }
         },
         "title": "NFCoreTemplateConfig",
         "type": "object"
      },
      "NFCoreYamlLintConfig": {
         "description": "schema for linting config in `.nf-core.yml` should cover:\n\n.. code-block:: yaml\n    files_unchanged:\n        - .github/workflows/branch.yml\n    modules_config: False\n    modules_config:\n            - fastqc\n    # merge_markers: False\n    merge_markers:\n            - docs/my_pdf.pdf\n    nextflow_config: False\n    nextflow_config:\n        - manifest.name\n        - config_defaults:\n            - params.annotation_db\n            - params.multiqc_comment_headers\n            - params.custom_table_headers\n    # multiqc_config: False\n    multiqc_config:\n        - report_section_order\n        - report_comment\n    files_exist:\n        - .github/CONTRIBUTING.md\n        - CITATIONS.md\n    template_strings: False\n    template_strings:\n            - docs/my_pdf.pdf\n    nfcore_components: False",
         "properties": {
            "files_unchanged": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Files Unchanged"
            },
            "modules_config": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Modules Config"
            },
            "merge_markers": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Merge Markers"
            },
            "nextflow_config": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "anyOf": [
                           {
                              "type": "string"
                           },
                           {
                              "additionalProperties": {
                                 "items": {
                                    "type": "string"
                                 },
                                 "type": "array"
                              },
                              "type": "object"
                           }
                        ]
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Nextflow Config"
            },
            "multiqc_config": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Multiqc Config"
            },
            "files_exist": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Files Exist"
            },
            "template_strings": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Template Strings"
            },
            "readme": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "items": {
                        "type": "string"
                     },
                     "type": "array"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Readme"
            },
            "nfcore_components": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Nfcore Components"
            },
            "actions_ci": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Actions Ci"
            },
            "actions_awstest": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Actions Awstest"
            },
            "actions_awsfulltest": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Actions Awsfulltest"
            },
            "pipeline_todos": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Pipeline Todos"
            },
            "pipeline_if_empty_null": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Pipeline If Empty Null"
            },
            "plugin_includes": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Plugin Includes"
            },
            "pipeline_name_conventions": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Pipeline Name Conventions"
            },
            "schema_lint": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Schema Lint"
            },
            "schema_params": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Schema Params"
            },
            "system_exit": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "System Exit"
            },
            "schema_description": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Schema Description"
            },
            "actions_schema_validation": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Actions Schema Validation"
            },
            "modules_json": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Modules Json"
            },
            "modules_structure": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Modules Structure"
            },
            "base_config": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Base Config"
            },
            "nfcore_yml": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Nfcore Yml"
            },
            "version_consistency": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Version Consistency"
            },
            "included_configs": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Included Configs"
            },
            "local_component_structure": {
               "anyOf": [
                  {
                     "type": "boolean"
                  },
                  {
                     "type": "null"
                  }
               ],
               "default": null,
               "title": "Local Component Structure"
            }
         },
         "title": "NFCoreYamlLintConfig",
         "type": "object"
      }
   }
}
```

</details></p>
* **Fields:**
  - [`bump_version (Dict[str, bool] | None)`](#nf_core.utils.NFCoreYamlConfig.bump_version)
  - [`lint (nf_core.utils.NFCoreYamlLintConfig | None)`](#nf_core.utils.NFCoreYamlConfig.lint)
  - [`nf_core_version (str | None)`](#nf_core.utils.NFCoreYamlConfig.nf_core_version)
  - [`org_path (str | None)`](#nf_core.utils.NFCoreYamlConfig.org_path)
  - [`repository_type (Literal['pipeline', 'modules'] | None)`](#nf_core.utils.NFCoreYamlConfig.repository_type)
  - [`template (nf_core.utils.NFCoreTemplateConfig | None)`](#nf_core.utils.NFCoreYamlConfig.template)
  - [`update (Dict[str, str | bool | Dict[str, str | Dict[str, str | bool]]] | None)`](#nf_core.utils.NFCoreYamlConfig.update)

#### _`field{:python}`_`bump_version{:python}`_: Dict\[str, bool] | None_`{:python}`_= None_

Disable bumping of the version for a module/subworkflow (when repository_type is modules). See <https://nf-co.re/docs/nf-core-tools/modules/bump-versions> for more information.

#### _`field{:python}`_`lint{:python}`_: [NFCoreYamlLintConfig](#nf_core.utils.NFCoreYamlLintConfig) | None_`{:python}`_= None_

Pipeline linting configuration, see <https://nf-co.re/docs/nf-core-tools/pipelines/lint#linting-config> for examples and documentation

#### _`field{:python}`_`nf_core_version{:python}`_: str | None_`{:python}`_= None_

Version of nf-core/tools used to create/update the pipeline

#### _`field{:python}`_`org_path{:python}`_: str | None_`{:python}`_= None_

Path to the organisation’s modules repository (used for modules repo_type only)

#### _`field{:python}`_`repository_type{:python}`_: Literal\['pipeline', 'modules'] | None_`{:python}`_= None_

Type of repository

#### _`field{:python}`_`template{:python}`_: [NFCoreTemplateConfig](#nf_core.utils.NFCoreTemplateConfig) | None_`{:python}`_= None_

Pipeline template configuration

#### _`field{:python}`_`update{:python}`_: Dict\[str, str | bool | Dict\[str, str | Dict\[str, str | bool]]] | None_`{:python}`_= None_

Disable updating specific modules/subworkflows (when repository_type is pipeline). See <https://nf-co.re/docs/nf-core-tools/modules/update> for more information.

#### `get(item: str, default: Any = None) → Any{:python}`

#### `model_dump(**kwargs) → Dict[str, Any]{:python}`

Usage docs: <https://docs.pydantic.dev/2.10/concepts/serialization/#modelmodel_dump>

Generate a dictionary representation of the model, optionally specifying which fields to include or exclude.

- **Parameters:**
  - **mode** – The mode in which to_python should run.
    If mode is ‘json’, the output will only contain JSON serializable types.
    If mode is ‘python’, the output may contain non-JSON-serializable Python objects.
  - **include** – A set of fields to include in the output.
  - **exclude** – A set of fields to exclude from the output.
  - **context** – Additional context to pass to the serializer.
  - **by_alias** – Whether to use the field’s alias in the dictionary key if defined.
  - **exclude_unset** – Whether to exclude fields that have not been explicitly set.
  - **exclude_defaults** – Whether to exclude fields that are set to their default value.
  - **exclude_none** – Whether to exclude fields that have a value of None.
  - **round_trip** – If True, dumped values should be valid as input for non-idempotent types such as Json\[T].
  - **warnings** – How to handle serialization errors. False/”none” ignores them, True/”warn” logs errors,
    “error” raises a \[PydanticSerializationError]\[pydantic_core.PydanticSerializationError].
  - **serialize_as_any** – Whether to serialize fields with duck-typing serialization behavior.
- **Returns:**
  A dictionary representation of the model.

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

### _`pydantic model{:python}`_`nf_core.utils.NFCoreYamlLintConfig{:python}`

Bases: `BaseModel`

schema for linting config in .nf-core.yml should cover:

<p><details  class="autodoc_pydantic_collapsable_json">
<summary>Show JSON schema</summary>
```json
{
   "title": "NFCoreYamlLintConfig",
   "description": "schema for linting config in `.nf-core.yml` should cover:\n\n.. code-block:: yaml\n    files_unchanged:\n        - .github/workflows/branch.yml\n    modules_config: False\n    modules_config:\n            - fastqc\n    # merge_markers: False\n    merge_markers:\n            - docs/my_pdf.pdf\n    nextflow_config: False\n    nextflow_config:\n        - manifest.name\n        - config_defaults:\n            - params.annotation_db\n            - params.multiqc_comment_headers\n            - params.custom_table_headers\n    # multiqc_config: False\n    multiqc_config:\n        - report_section_order\n        - report_comment\n    files_exist:\n        - .github/CONTRIBUTING.md\n        - CITATIONS.md\n    template_strings: False\n    template_strings:\n            - docs/my_pdf.pdf\n    nfcore_components: False",
   "type": "object",
   "properties": {
      "files_unchanged": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Files Unchanged"
      },
      "modules_config": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Modules Config"
      },
      "merge_markers": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Merge Markers"
      },
      "nextflow_config": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "anyOf": [
                     {
                        "type": "string"
                     },
                     {
                        "additionalProperties": {
                           "items": {
                              "type": "string"
                           },
                           "type": "array"
                        },
                        "type": "object"
                     }
                  ]
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Nextflow Config"
      },
      "multiqc_config": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Multiqc Config"
      },
      "files_exist": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Files Exist"
      },
      "template_strings": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Template Strings"
      },
      "readme": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Readme"
      },
      "nfcore_components": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Nfcore Components"
      },
      "actions_ci": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Actions Ci"
      },
      "actions_awstest": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Actions Awstest"
      },
      "actions_awsfulltest": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Actions Awsfulltest"
      },
      "pipeline_todos": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Pipeline Todos"
      },
      "pipeline_if_empty_null": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Pipeline If Empty Null"
      },
      "plugin_includes": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Plugin Includes"
      },
      "pipeline_name_conventions": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Pipeline Name Conventions"
      },
      "schema_lint": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Schema Lint"
      },
      "schema_params": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Schema Params"
      },
      "system_exit": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "System Exit"
      },
      "schema_description": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Schema Description"
      },
      "actions_schema_validation": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Actions Schema Validation"
      },
      "modules_json": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Modules Json"
      },
      "modules_structure": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Modules Structure"
      },
      "base_config": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Base Config"
      },
      "nfcore_yml": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Nfcore Yml"
      },
      "version_consistency": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Version Consistency"
      },
      "included_configs": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Included Configs"
      },
      "local_component_structure": {
         "anyOf": [
            {
               "type": "boolean"
            },
            {
               "type": "null"
            }
         ],
         "default": null,
         "title": "Local Component Structure"
      }
   }
}
```

</details></p>
* **Fields:**
  - [`actions_awsfulltest (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.actions_awsfulltest)
  - [`actions_awstest (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.actions_awstest)
  - [`actions_ci (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.actions_ci)
  - [`actions_schema_validation (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.actions_schema_validation)
  - [`base_config (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.base_config)
  - [`files_exist (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.files_exist)
  - [`files_unchanged (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.files_unchanged)
  - [`included_configs (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.included_configs)
  - [`local_component_structure (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.local_component_structure)
  - [`merge_markers (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.merge_markers)
  - [`modules_config (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.modules_config)
  - [`modules_json (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.modules_json)
  - [`modules_structure (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.modules_structure)
  - [`multiqc_config (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.multiqc_config)
  - [`nextflow_config (bool | List[str | Dict[str, List[str]]] | None)`](#nf_core.utils.NFCoreYamlLintConfig.nextflow_config)
  - [`nfcore_components (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.nfcore_components)
  - [`nfcore_yml (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.nfcore_yml)
  - [`pipeline_if_empty_null (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.pipeline_if_empty_null)
  - [`pipeline_name_conventions (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.pipeline_name_conventions)
  - [`pipeline_todos (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.pipeline_todos)
  - [`plugin_includes (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.plugin_includes)
  - [`readme (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.readme)
  - [`schema_description (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.schema_description)
  - [`schema_lint (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.schema_lint)
  - [`schema_params (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.schema_params)
  - [`system_exit (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.system_exit)
  - [`template_strings (bool | List[str] | None)`](#nf_core.utils.NFCoreYamlLintConfig.template_strings)
  - [`version_consistency (bool | None)`](#nf_core.utils.NFCoreYamlLintConfig.version_consistency)

#### _`field{:python}`_`actions_awsfulltest{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to run full tests on AWS

#### _`field{:python}`_`actions_awstest{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to run tests on AWS

#### _`field{:python}`_`actions_ci{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to use GitHub Actions CI

#### _`field{:python}`_`actions_schema_validation{:python}`_: bool | None_`{:python}`_= None_

Lint GitHub Action workflow files with schema

#### _`field{:python}`_`base_config{:python}`_: bool | None_`{:python}`_= None_

Lint base.config file

#### _`field{:python}`_`files_exist{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of files that can not exist

#### _`field{:python}`_`files_unchanged{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of files that should not be changed

#### _`field{:python}`_`included_configs{:python}`_: bool | None_`{:python}`_= None_

Lint for included configs

#### _`field{:python}`_`local_component_structure{:python}`_: bool | None_`{:python}`_= None_

Lint local components use correct structure mirroring remote

#### _`field{:python}`_`merge_markers{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of files that should not contain merge markers

#### _`field{:python}`_`modules_config{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of modules that should not be changed

#### _`field{:python}`_`modules_json{:python}`_: bool | None_`{:python}`_= None_

Lint modules.json file

#### _`field{:python}`_`modules_structure{:python}`_: bool | None_`{:python}`_= None_

Lint modules structure

#### _`field{:python}`_`multiqc_config{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of MultiQC config options that be changed

#### _`field{:python}`_`nextflow_config{:python}`_: bool | List\[str | Dict\[str, List\[str]]] | None_`{:python}`_= None_

List of Nextflow config files that should not be changed

#### _`field{:python}`_`nfcore_components{:python}`_: bool | None_`{:python}`_= None_

Lint all required files to use nf-core modules and subworkflows

#### _`field{:python}`_`nfcore_yml{:python}`_: bool | None_`{:python}`_= None_

Lint nf-core.yml

#### _`field{:python}`_`pipeline_if_empty_null{:python}`_: bool | None_`{:python}`_= None_

Lint for ifEmpty(null) statements

#### _`field{:python}`_`pipeline_name_conventions{:python}`_: bool | None_`{:python}`_= None_

Lint for pipeline name conventions

#### _`field{:python}`_`pipeline_todos{:python}`_: bool | None_`{:python}`_= None_

Lint for TODOs statements

#### _`field{:python}`_`plugin_includes{:python}`_: bool | None_`{:python}`_= None_

Lint for nextflow plugin

#### _`field{:python}`_`readme{:python}`_: bool | List\[str] | None_`{:python}`_= None_

Lint the README.md file

#### _`field{:python}`_`schema_description{:python}`_: bool | None_`{:python}`_= None_

Check that every parameter in the schema has a description.

#### _`field{:python}`_`schema_lint{:python}`_: bool | None_`{:python}`_= None_

Lint nextflow_schema.json file

#### _`field{:python}`_`schema_params{:python}`_: bool | None_`{:python}`_= None_

Lint schema for all params

#### _`field{:python}`_`system_exit{:python}`_: bool | None_`{:python}`_= None_

Lint for System.exit calls in groovy/nextflow code

#### _`field{:python}`_`template_strings{:python}`_: bool | List\[str] | None_`{:python}`_= None_

List of files that can contain template strings

#### _`field{:python}`_`version_consistency{:python}`_: bool | None_`{:python}`_= None_

Lint for version consistency

#### `get(item: str, default: Any = None) → Any{:python}`

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

### _`class{:python}`_`nf_core.utils.Pipeline(wf_path: Path){:python}`

Bases: `object`

Object to hold information about a local pipeline.

- **Parameters:**
  **path** (_str_) – The path to the nf-core pipeline directory.

#### `conda_config{:python}`

The parsed conda configuration file content (`environment.yml`).

- **Type:**
  dict

#### `conda_package_info{:python}`

The conda package(s) information, based on the API requests to Anaconda cloud.

- **Type:**
  dict

#### `nf_config{:python}`

The Nextflow pipeline configuration file content.

- **Type:**
  dict

#### `files{:python}`

A list of files found during the linting process.

- **Type:**
  list

#### `git_sha{:python}`

The git sha for the repo commit / current GitHub pull-request ($GITHUB_PR_COMMIT)

- **Type:**
  str

#### `minNextflowVersion{:python}`

The minimum required Nextflow version to run the pipeline.

- **Type:**
  str

#### `wf_path{:python}`

Path to the pipeline directory.

- **Type:**
  str

#### `pipeline_name{:python}`

The pipeline name, without the nf-core tag, for example hlatyping.

- **Type:**
  str

#### `schema_obj{:python}`

A `PipelineSchema` object

- **Type:**
  obj

#### `_fp(fn: str | Path) → Path{:python}`

Convenience function to get full path to a file in the pipeline

#### `_load() → bool{:python}`

Run core load functions

#### `_load_conda_environment() → bool{:python}`

Try to load the pipeline environment.yml file, if it exists

#### `list_files() → List[Path]{:python}`

Get a list of all files in the pipeline

#### `load_pipeline_config() → bool{:python}`

Get the nextflow config for this pipeline

Once loaded, set a few convenience reference class attributes

### _`class{:python}`_`nf_core.utils.SingularityCacheFilePathValidator{:python}`

Bases: `Validator`

Validator for file path specified as –singularity-cache-index argument in nf-core pipelines download

#### `_abc_impl{:python}`_= <\_abc.\_abc_data object>_

#### `validate(value){:python}`

Validate the input.
If invalid, this should raise a `ValidationError`.

- **Parameters:**
  **document** – `Document` instance.

### `nf_core.utils.anaconda_package(dep, dep_channels=None){:python}`

Query conda package information.

Sends a HTTP GET request to the Anaconda remote API.

- **Parameters:**
  - **dep** (_str_) – A conda package name.
  - **dep_channels** (_list_) – list of conda channels to use
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(\*\***404\***\*)** –

### `nf_core.utils.check_if_outdated(current_version=None, remote_version=None, source_url='https://nf-co.re/tools_version'){:python}`

Check if the current version of nf-core is outdated

### `nf_core.utils.custom_yaml_dumper(){:python}`

Overwrite default PyYAML output to make Prettier YAML linting happy

### `nf_core.utils.determine_base_dir(directory: Path | str = '.') → Path{:python}`

### `nf_core.utils.fetch_remote_version(source_url){:python}`

### `nf_core.utils.fetch_wf_config(wf_path: Path, cache_config: bool = True) → dict{:python}`

Uses Nextflow to retrieve the the configuration variables
from a Nextflow workflow.

- **Parameters:**
  - **wf_path** (_str_) – Nextflow workflow file system path.
  - **cache_config** (_bool_) – cache configuration or not (def. True)
- **Returns:**
  Workflow configuration settings.
- **Return type:**
  dict

### `nf_core.utils.file_md5(fname){:python}`

Calculates the md5sum for a file on the disk.

- **Parameters:**
  **fname** (_str_) – Path to a local file.

### `nf_core.utils.get_biocontainer_tag(package, version){:python}`

Given a bioconda package and version, looks for Docker and Singularity containers
using the biocontaineres API, e.g.:
<https://api.biocontainers.pro/ga4gh/trs/v2/tools>/{tool}/versions/{tool}-{version}
Returns the most recent container versions by default.
:param package: A bioconda package name.
:type package: str
:param version: Version of the bioconda package
:type version: str

- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** **or** **gives an unexpected status code** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** **(\*\***404\***\*)** –

### `nf_core.utils.get_first_available_path(directory: Path | str, paths: List[str]) → Path | None{:python}`

### `nf_core.utils.get_repo_commit(pipeline, commit_id){:python}`

Check if the repo contains the requested commit_id, and expand it to long form if necessary.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **commit_id** – The requested commit ID (SHA). It can be in standard long/short form, or any length.
- **Returns:**
  String or None
- **Return type:**
  commit_id

### `nf_core.utils.get_repo_releases_branches(pipeline, wfs){:python}`

Fetches details of a nf-core workflow to download.

- **Parameters:**
  - **pipeline** (_str_) – GitHub repo username/repo
  - **wfs** – A nf_core.pipelines.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  Array of releases, Array of branches
- **Return type:**
  wf_releases, wf_branches (tuple)
- **Raises:**
  **LockupError**\*\*,\*\* **if the pipeline can not be found.** –

### `nf_core.utils.get_wf_files(wf_path: Path){:python}`

Return a list of all files in a directory (ignores .gitigore files)

### `nf_core.utils.is_file_binary(path){:python}`

Check file path to see if it is a binary file

### `nf_core.utils.is_pipeline_directory(wf_path){:python}`

Checks if the specified directory have the minimum required files
(‘main.nf’, ‘nextflow.config’) for a pipeline directory

- **Parameters:**
  **wf_path** (_str_) – The directory to be inspected
- **Raises:**
  **UserWarning** – If one of the files are missing

### `nf_core.utils.is_relative_to(path1, path2){:python}`

Checks if a path is relative to another.

Should mimic Path.is_relative_to which not available in Python < 3.9

path1 (Path | str): The path that could be a subpath
path2 (Path | str): The path the could be the superpath

### `nf_core.utils.load_tools_config(directory: str | Path = '.') → Tuple[Path | None,{:python}`[`NFCoreYamlConfig{:python}`](#nf_core.utils.NFCoreYamlConfig)`| None]{:python}`

Parse the nf-core.yml configuration file

Look for a file called either .nf-core.yml or .nf-core.yaml

Also looks for the deprecated file .nf-core-lint.yml/yaml and issues
a warning that this file will be deprecated in the future

Returns the loaded config dict or False, if the file couldn’t be loaded

### `nf_core.utils.nested_delitem(d, keys){:python}`

Deletes a key from a nested dictionary

- **Parameters:**
  - **d** (_dict_) – the nested dictionary to traverse
  - **keys** (_list_ \*\[\*_Any_ _]_) – A list of keys to iteratively traverse, deleting the final one

### `nf_core.utils.nested_setitem(d, keys, value){:python}`

Sets the value in a nested dict using a list of keys to traverse

- **Parameters:**
  - **d** (_dict_) – the nested dictionary to traverse
  - **keys** (_list_ \*\[\*_Any_ _]_) – A list of keys to iteratively traverse
  - **value** (_Any_) – The value to be set for the last key in the chain

### `nf_core.utils.parse_anaconda_licence(anaconda_response, version=None){:python}`

Given a response from the anaconda API using anaconda_package, parse the software licences.

Returns: Set of licence types

### `nf_core.utils.pip_package(dep){:python}`

Query PyPI package information.

Sends a HTTP GET request to the PyPI remote API.

- **Parameters:**
  **dep** (_str_) – A PyPI package name.
- **Raises:**
  - **A LookupError**\*\*,\*\* **if the connection fails** **or** **times out** –
  - **A ValueError**\*\*,\*\* **if the package name can not be found** –

### `nf_core.utils.plural_es(list_or_int){:python}`

Return a ‘es’ if the input is not one or has not the length of one.

### `nf_core.utils.plural_s(list_or_int){:python}`

Return an s if the input is not one or has not the length of one.

### `nf_core.utils.plural_y(list_or_int){:python}`

Return ‘ies’ if the input is not one or has not the length of one, else ‘y’.

### `nf_core.utils.poll_nfcore_web_api(api_url: str, post_data: Dict | None = None) → Dict{:python}`

Poll the nf-core website API

Takes argument api_url for URL

Expects API response to be valid JSON and contain a top-level ‘status’ key.

### `nf_core.utils.prompt_pipeline_release_branch(wf_releases: List[Dict[str, Any]], wf_branches: Dict[str, Any], multiple: bool = False) → Tuple[Any, List[str]]{:python}`

Prompt for pipeline release / branch

- **Parameters:**
  - **wf_releases** (_array_) – Array of repo releases as returned by the GitHub API
  - **wf_branches** (_array_) – Array of repo branches, as returned by the GitHub API
  - **multiple** (_bool_) – Allow selection of multiple releases & branches (for Seqera Platform)
- **Returns:**
  Selected release / branch or False if no releases / branches available
- **Return type:**
  choice (questionary.Choice or bool)

### `nf_core.utils.prompt_remote_pipeline_name(wfs){:python}`

Prompt for the pipeline name with questionary

- **Parameters:**
  **wfs** – A nf_core.pipelines.list.Workflows() object, where get_remote_workflows() has been called.
- **Returns:**
  GitHub repo - username/repo
- **Return type:**
  pipeline (str)
- **Raises:**
  **AssertionError**\*\*,\*\* **if pipeline cannot be found** –

### `nf_core.utils.rich_force_colors(){:python}`

Check if any environment variables are set to force Rich to use coloured output

### `nf_core.utils.run_cmd(executable: str, cmd: str) → Tuple[bytes, bytes] | None{:python}`

Run a specified command and capture the output. Handle errors nicely.

### `nf_core.utils.set_wd(path: Path) → Generator[None, None, None]{:python}`

Sets the working directory for this context.

- **Parameters:**
  **path** (_Path_) – Path to the working directory to be used inside this context.

### `nf_core.utils.setup_nfcore_cachedir(cache_fn: str | Path) → Path{:python}`

Sets up local caching for caching files between sessions.

### `nf_core.utils.setup_nfcore_dir() → bool{:python}`

Creates a directory for files that need to be kept between sessions

Currently only used for keeping local copies of modules repos

### `nf_core.utils.setup_requests_cachedir() → Dict[str, Path | timedelta | str]{:python}`

Sets up local caching for faster remote HTTP requests.

Caching directory will be set up in the user’s home directory under
a .config/nf-core/cache\_\* subdir.

Uses requests_cache monkey patching.
Also returns the config dict so that we can use the same setup with a Session.

### `nf_core.utils.sort_dictionary(d: Dict) → Dict{:python}`

Sorts a nested dictionary recursively

### `nf_core.utils.strip_ansi_codes(string, replace_with=''){:python}`

Strip ANSI colouring codes from a string to return plain text.

From Stack Overflow: <https://stackoverflow.com/a/14693789/713980>

### `nf_core.utils.validate_file_md5(file_name, expected_md5hex){:python}`

Validates the md5 checksum of a file on disk.

- **Parameters:**
  - **file_name** (_str_) – Path to a local file.
  - **expected** (_str_) – The expected md5sum.
- **Raises:**
  **IOError**\*\*,\*\* **if the md5sum does not match the remote sum.** –

### `nf_core.utils.wait_cli_function(poll_func: Callable[[], bool], refresh_per_second: int = 20) → None{:python}`

Display a command-line spinner while calling a function repeatedly.

Keep waiting until that function returns True

- **Parameters:**
  - **poll_func** (_function_) – Function to call
  - **refresh_per_second** (_int_) – Refresh this many times per second. Default: 20.
- **Returns:**
  None. Just sits in an infinite loop until the function returns True.
