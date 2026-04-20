---
title: Reviewing nf-core configs
subtitle: Review configs pull requests
shortTitle: Reviewing configs
markdownPlugin: checklist
---

nf-core/configs reviews ensure that [configs](https://nf-co.re/configs/).
When you review a nf-core/configs pull request, you examine a new config submission or proposed changes to an existing component and provide constructive feedback before maintainers merge them into the nf-core repository.

:::note{title="Component specifications"}
For information about best practices for nf-core configs, see [configs specifications](../../../specifications/configs/overview#configs-specifications).
:::

## General

- [ ] Permission by infrastructure administrators confirmed

For institutional configs:

- [ ] Each config has a `*.conf` Nextflow configuration file under `conf/`
- [ ] Each file has a `*md` documentation file under `docs/`
- [ ] Does not include `withName` or `withLabel` specifications
- [ ] Includes `resourceLimits` scope

For pipeline-specific institutional configs

- [ ] Each config targets a single pipeline
- [ ] Each config has a `*.conf` Nextflow configuration file under `conf/pipeline/<pipele_name>/`
- [ ] Each file has a `*md` documentation file under `docs/pipeline/<pipeline name>/`

## Naming

- [ ] The name is all lowercase
- [ ] The name includes only letters or numbers
- [ ] The name includes only underscore
- [ ] The name uses no other characters or punctuation
- [ ] If multiple clusters, uses an institutional suffix

## Config contents

- [ ] Each config adheres to the Nextflow strict syntax
- [ ] Each config includes parameters:
  - [ ] `config_profile_description`
  - [ ] `config_profile_contact`
  - [ ] `config_profile_url`
- [ ] Custom parameters:
  - [ ] Are documented
  - [ ] Are added to the `ignoreParams` option of the `validation` scope

## Documentation

- [ ] Contains information on how to use the config in the context of the instruction
- [ ] (Optional) documentation includes additional recommendations (e.g. scratch, temporary, or environmental variables)
- [ ] Required parameter or environmental variables described
- [ ] All commands are copy-paste-able examples

## Config infrastructure

- [ ] Entry added to `nfcore_custom.config`
- [ ] Entry added to `.github/workflows/main.yaml`
