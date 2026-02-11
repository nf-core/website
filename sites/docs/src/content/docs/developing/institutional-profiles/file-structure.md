---
title: File structure requirements
subtitle: Which files to create or modify for your institutional profile
shortTitle: File structure
---

Institutional profiles require specific files in the [nf-core/configs](https://github.com/nf-core/configs) repository.
You will create new files for your cluster configuration and documentation, and modify existing files to register your profile.

## Files to create

Create these new files for your institutional profile:

### Configuration file

**Path**: `conf/<cluster_name>.config`

This is the main configuration file that defines your cluster's settings.
It contains all the Nextflow configuration scopes (parameters, process settings, executor options, container configurations) that apply to your cluster.

Use lowercase with underscores for the filename (for example, `big_university.config`, `research_cluster.config`).

See [Configuration file components](configuration.md) for detailed information about what to include.

### Documentation file

**Path**: `docs/<cluster_name>.md`

This markdown file documents your profile's purpose, usage instructions, and any special requirements.
It helps users at your institution understand how to use the profile and troubleshoot common issues.

The filename should match your configuration file name (for example, `big_university.md` for `big_university.config`).

See [Documentation requirements](documentation.md) for detailed information about what to include.

## Files to modify

Update these existing files in the nf-core/configs repository to register your profile:

### nfcore_custom.config

**Path**: `nfcore_custom.config`

Add a profile entry that references your configuration file.
This makes your profile available when users specify `-profile <your_institution>`.

Add your profile within the `profiles` scope in alphabetical order:

```groovy
profiles {
  // ... other profiles ...
  big_university {
    includeConfig "${params.custom_config_base}/conf/big_university.config"
  }
  // ... other profiles ...
}
```

The profile name in this file determines how users access your profile with the `-profile` flag.

### README.md

**Path**: `README.md`

Add your institution to the documentation table in the main README.
This provides visibility for users browsing the repository.

Add a row to the institution table in alphabetical order:

```markdown
| big_university | Big University HPC | [big_university.md](docs/big_university.md) |
```

Include the profile name, institution name, and link to your documentation file.

### GitHub Actions workflow

**Path**: `.github/workflows/main.yml`

Add your profile to the automated testing matrix.
This ensures your profile passes continuous integration checks when template changes occur.

Locate the `profile` list in the test matrix and add your profile name in alphabetical order:

```yaml
strategy:
  matrix:
    profile:
      - abims
      - awsbatch
      - big_university # Add your profile here
      - bigpurple
      # ... other profiles ...
```

:::note
The GitHub Actions workflow file may have a different structure depending on when you create your profile.
Follow the existing pattern in the file and add your profile to the appropriate location in the test matrix.
:::

## Pipeline-specific profiles

If you need a pipeline-specific profile (not common), create an additional directory and configuration file:

**Path**: `conf/pipeline/<pipeline_name>/<cluster_name>.config`

For example, to create a profile specific to the rnaseq pipeline for big_university:

```
conf/pipeline/rnaseq/big_university.config
```

Then reference it in `nfcore_custom.config`:

```groovy
profiles {
  big_university {
    includeConfig "${params.custom_config_base}/conf/big_university.config"
    includeConfig "${params.custom_config_base}/conf/pipeline/rnaseq/big_university.config"
  }
}
```

The pipeline-specific configuration loads after the global profile, so it can override or extend global settings.

## Next steps

Once you understand which files to create and modify, you can begin writing your configuration file.

See [Configuration file components](configuration.md) to learn how to structure your cluster configuration.
