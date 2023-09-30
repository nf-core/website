---
title: 'Tutorial: Using nf-core configs outside nf-core'
subtitle: Guidance on how to use nf-core's centralised institutional configs in your own workflows.
---

When you've been working with nf-core pipelines for a while, you often will get very used to many of the convenience functionality offered by the wider nf-core infrastrucutre.

One such thing is the centralised nf-core/configs repository of pre-configured configuration files that allow nf-core pipelines to run optimally on institutional clusters via the `-profile` parameter, e.g. `-profile uppmax`. A list of existing institutional profiles can be seen on the [nf-core website](https://nf-co.re/configs).

One great thing about nf-core/configs is that they aren't just restricted to nf-core pipelines, they can also be used in [fully fledged 'unofficial' nf-core pipelines](unofficial_pipelines.md) but also in your own custom mini-scripts and pipelines!

Here we will describe the four things you will need to do in your custom script or pipeline to use nf-core institutional configs.

1. Set default basic resources for all processes, e.g. in a `conf/base.config`

At the bare minimum, the file should contain:

```groovy
process {
    cpus   = { check_max( 1    * task.attempt, 'cpus'   ) }
    memory = { check_max( 7.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h  * task.attempt, 'time'   ) }
}
```

(or other sensible default values)

For a more sophisticated `base.config`, see the full [nf-core template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config)

2. Load the base config in the `nextflow.config`

The following should be placed outside all scopes:

```groovy
includeConfig 'conf/base.config'
```

3. Load nf-core's institutional profile repository in the pipeline's `nextflow.config`

```groovy
// Load nf-core custom profiles from different Institutions
try {
    includeConfig "https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config"
} catch (Exception e) {
    System.err.println("WARNING: Could not load nf-core/config profiles: https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config")
}
```

Note nf-core pipelines use a [parameter for the URL](https://github.com/nf-core/tools/blob/0912990a63ef29e44e07cc2ba6ab81113684e0ae/nf_core/pipeline-template/nextflow.config#L67-L72) to allow debugging.

4. Finally Add the the [`check_max`](https://github.com/nf-core/tools/blob/0912990a63ef29e44e07cc2ba6ab81113684e0ae/nf_core/pipeline-template/nextflow.config#L233-L264) function to `nextflow.config`

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

With this, you should be able to run `nextflow run mainf.nf -profile <your_institutional_profile>`, and your custom script/pipeline should integrate nicely with your cluster!
