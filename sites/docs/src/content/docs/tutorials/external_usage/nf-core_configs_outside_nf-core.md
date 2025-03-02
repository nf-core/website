---
title: "Using nf-core configs outside nf-core"
subtitle: Guidance on how to use nf-core's centralised institutional configs in your own workflows.
---

When you've been working with nf-core pipelines for a while, you often will get very used to many of the convenience functionality offered by the wider nf-core infrastructure.

One such thing is the centralised nf-core/configs repository of pre-configured configuration files that allow nf-core pipelines to run optimally on institutional clusters via the `-profile` parameter, e.g. `-profile uppmax`. A list of existing institutional profiles can be seen on the [nf-core website](https://nf-co.re/configs).

:::tip
If you want to write your own institutional profile, see the [guide on how to write a new institutional profile](/docs/tutorials/use_nf-core_pipelines/writing_institutional_profiles).
:::

One great thing about nf-core/configs is that they aren't just restricted to nf-core pipelines, they can also be used in [fully fledged 'unofficial' nf-core pipelines](/docs/guidelines/external_use) but also in your own custom mini-scripts and pipelines!

Here we will describe the steps you will need to perform in your custom script or pipeline to use nf-core institutional configs.

1. In a `conf/base.config` file set default basic resources for all processes. Can be named something else.

   At the bare minimum, the file should contain:

   ```groovy
   process {
       cpus   = 1
       memory = 7.GB
       time   = 4.h
   }
   ```

   (or other sensible default values)

   :::note{collapse title="Note on older nf-core template/Nextflow versions"}

   When running pipelines generated with the nf-core template before version v3.0.0 or with Nextflow before version 24.04.0 you may need to use an older syntax for setting resources. The following closures prevent resources from exceeding a maximum limit.

   Set in `conf/base.config` default values.

   ```groovy
   process {
       cpus   = { check_max( 1    * task.attempt, 'cpus'   ) }
       memory = { check_max( 7.GB * task.attempt, 'memory' ) }
       time   = { check_max( 4.h  * task.attempt, 'time'   ) }
   }
   ```

   Then add the [`check_max`](https://github.com/nf-core/tools/blob/0912990a63ef29e44e07cc2ba6ab81113684e0ae/nf_core/pipeline-template/nextflow.config#L233-L264) function to `nextflow.config`

   The following should be placed outside all scopes, typically at the bottom of the file:

   ```groovy
   // Function to ensure that resource requirements don't go beyond
   // a maximum limit
   def check_max(obj, type) {
     if (type == 'memory') {
         try {
             if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                 return params.max_memory as nextflow.util.MemoryUnit
             else
                 return obj
         } catch (all) {
             println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
             return obj
         }
     } else if (type == 'time') {
         try {
             if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                 return params.max_time as nextflow.util.Duration
             else
                 return obj
         } catch (all) {
             println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
             return obj
         }
     } else if (type == 'cpus') {
         try {
             return Math.min( obj, params.max_cpus as int )
         } catch (all) {
             println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
             return obj
         }
     }
   }
   ```

   :::

   For a more sophisticated `base.config`, see the full [nf-core template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config)

2. In a top level `nextflow.config`, specify two `params` to specify the URL where to look for nf-core/configs from.

   ```groovy
   params {
     custom_config_version      = 'master'
     custom_config_base         = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
   }
   ```

   The two parameters together make it easy to test from forks and development profiles on specific branches.

3. In the top level `nextflow.config`, load the newly made base config.

   The following should be placed outside all scopes:

   ```groovy
   includeConfig 'conf/base.config'
   ```

4. In the top level `nextflow.config`, load nf-core's institutional profile repository based on the new `params`. This should be placed _after_ the `conf/base.config` include.

   ```groovy
   includeConfig !System.getenv('NXF_OFFLINE') && params.custom_config_base ? "${params.custom_config_base}/nfcore_custom.config" : "/dev/null"
   ```

   Note that instead of the parameters, you can just directly specify a URL, e.g.

   ```groovy
   includeConfig !System.getenv('NXF_OFFLINE') ? 'https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config' : "/dev/null"
   ```

   :::note{collapse title="Note on older nf-core template/Nextflow versions"}
   If you wish to use pipelines generated with the nf-core template before `v3.0.0` , and/or when running with Nextflow versions earlier than 24.04.0 you may need to use an older syntax for loading the configs:

   ```groovy
   // Load nf-core custom profiles from different Institutions
   try {
       includeConfig "https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config"
   } catch (Exception e) {
       System.err.println("WARNING: Could not load nf-core/config profiles: https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config")
   }
   ```

   :::

With this, you should be able to run `nextflow run mainf.nf -profile <your_institutional_profile>`, and your custom script/pipeline should integrate nicely with your cluster!
