---
title: Using nf-core configs outside nf-core
subtitle: Learn how to integrate institutional profiles into custom workflows
shortTitle: Profiles outside nf-core
weight: 7
---

Institutional profiles from the nf-core/configs repository can be integrated into custom scripts and workflows, extending their utility beyond nf-core pipelines. These centralized configurations provide pre-configured cluster settings without duplicating configuration management across multiple projects.

## Integration steps

To integrate nf-core configs outside nf-core:

1. Create a `conf/base.config` file with default process resources:

   ```groovy
   process {
       cpus   = 1
       memory = 7.GB
       time   = 4.h
   }
   ```

   :::note
   Older nf-core templates (pre-v3.0.0) and Nextflow versions before 24.04.0 require alternative syntax using the `check_max()` function to enforce resource limits.
   :::

1. In your top-level `nextflow.config`, define parameters for accessing the institutional configs repository:

   ```groovy
   params {
     custom_config_version      = 'master'
     custom_config_base         = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
   }
   ```

1. Include the base configuration file in your `nextflow.config`:

   ```groovy
   includeConfig 'conf/base.config'
   ```

1. Load the nf-core institutional repository:

   ```groovy
   includeConfig !System.getenv('NXF_OFFLINE') && params.custom_config_base ? "${params.custom_config_base}/nfcore_custom.config" : "/dev/null"
   ```

   :::note
   This conditional statement:

   - Checks if Nextflow is running in offline mode
   - Verifies the `custom_config_base` parameter exists
   - Loads the nf-core configs if both conditions are met
   - Falls back to `/dev/null` if offline or the parameter is not set
     :::

## Alternative approach

For simpler setups, use a direct URL approach without conditional logic:

```groovy
includeConfig 'https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config'
```

:::warning
This approach does not support offline mode or testing against development branches.
:::

## Usage

Once configured, execute your custom workflow using an institutional profile:

```bash
nextflow run main.nf -profile <institutional_profile>
```

Replace `<institutional_profile>` with your institution's profile name (for example, `biohpc_gen`, `uppmax`, `czbiohub_aws`).

You can also combine multiple profiles:

```bash
nextflow run main.nf -profile <institutional_profile>,docker
```

## Troubleshooting

### Profile not found

**Problem:** You receive a "Profile not found" error when running your workflow

**Solution:** Verify your profile configuration and network access:

- Check that the profile name matches an existing profile in nf-core/configs
- Verify the `custom_config_base` parameter points to a valid URL
- Ensure you have network access to the nf-core/configs repository
- Test the URL directly in a browser to confirm it is accessible

### Offline mode

**Problem:** Institutional profiles do not load when running in offline mode

**Solution:** Configure your workflow to use local profile files:

1. Download the required profile files from `nf-core/configs` to your local system
2. Update the `includeConfig` statement in `nextflow.config` to point to the local file path
3. Remove the conditional check for `NXF_OFFLINE` to always load from local files
4. Ensure the local files remain synchronized with updates to nf-core/configs if needed

### Resource limits not enforced

**Problem:** Resource limits are not being enforced correctly in your workflow

**Solution:** Check your configuration and Nextflow version:

- Verify your Nextflow version (version 24.04.0 or later recommended for `resourceLimits` support)
- Confirm the institutional profile includes `resourceLimits` in the process scope
- Check that your base configuration does not override profile settings (load base config before institutional profiles)
- Review the order of `includeConfig` statements in `nextflow.config`
