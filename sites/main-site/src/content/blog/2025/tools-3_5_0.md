---
title: "nf-core/tools - 3.5.0"
subtitle: "Hot topic release"
headerImage: "https://images.unsplash.com/photo-1696192228333-e108d82a7829"
headerImageAlt: "Photo of a motorbike full with plastic bags filled with water and gold fish."
pubDate: 2025-11-19T12:00:00+02:00
authors:
  - "mashehu"
  - "mirpedrol"
label:
  - "tools"
---

This release is small in scope, but introduces major changes to the syntax in nf-core/modules.
These changes are part of the comming Nextflow syntax changes, which will be delivered to pipelines through template updates gradually.
For more information on the Nextflow syntax adoption in nf-core pipelines, you can read [the blogpost](https://nf-co.re/blog/2025/nextflow_syntax_nf-core_roadmap) detailing the roadmap.

## Topic channels for version handling in modules

With a heroic push during the hackathon led by [@nvnieuwk](https://github.com/nvnieuwk), we switched to [topic channels](https://www.nextflow.io/docs/latest/reference/channel.html#topic) for version handling in nf-core/modules.
This means we don't write the tool versions to a `versions.yml` file anymore, but instead use simpler channel logic to broadcast and collect the software versions of the tools used in modules.

The main change happens in the `main.nf` files:

```groovy title="main.nf"
output:
tuple val(meta), path("*.html"), emit: html
tuple val(meta), path("*.zip") , emit: zip
path  "versions.yml"           , emit: versions // [!code --]
tuple val("${task.process}"), val('fastqc'), eval('fastqc --version | sed "/FastQC v/!d; s/.*v//"'), emit: versions_fastqc, topic: versions // [!code ++]
︙
cat <<-END_VERSIONS > versions.yml // [!code --]
"${task.process}": // [!code --]
    fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' ) // [!code --]
END_VERSIONS // [!code --]
```

We updated the modules template, so if you run `nf-core modules create{:bash}` you will get the new syntax.
For now nf-core linting accepts both ways to collect versions, but gives a warning if you use the old syntax.
In your pipeline you can mix and match modules with the new and old syntax.

:::note
Note that the channel called `versions` is now renamed to `versions_<TOOL NAME>`, so you will have to adapt the channel names if you want to continue using the old logic for collection versions.
:::
We will slowly migrate all nf-core modules to the new syntax.

For more information on how to migrate and handle these new modules, please refer to the [migration guide](/docs/tutorials/migrate_to_topics/update_modules).
For more information about topic channels in general, we have a dedicated [blogpost about it](/blog/2025/version_topics).

We outlined future adoptions of new Nextflow syntax elements in the nf-core infrastructure talk during the Nextflow summit 20205 and in an upcoming blog post.

## AWS bug fix

With nf-core/tools v3.4.0, we introduced a bug in the `nextflow.config` for pipelines using iGenomes, leading to permission errors when querying AWS S3 buckets.
This has been fixed in this release.

## Changes on sync PRs

When a new release of nf-core/tools is out, all nf-core pipelines receive a sync PR with the updates in the template.
Before, these sync PRs were closed if they were not merged when a new one was opened.
Now, we will keep all PRs open. This will make it easier to do incremental template updates, which will help bringing older pipelines up to date.

As a new detail, the sync PRs will also include a link to the tools release blogpost,
to make it easier for pipeline maintainers to access documentation about the updates.

To learn more about these decisions, you can read the above mentioned [blogpost about the Nextflow syntax roadmap for nf-core](https://nf-co.re/blog/2025/nextflow_syntax_nf-core_roadmap).

## Patch release 3.5.1

Before we actually started the automated template sync, we found a small bug in the updated version of `nf-core pipelines sync`, which was fixed in 3.5.1.
These are the only changes in this patch release.

## Changelog

You can find the complete changelog and technical details [on GitHub](https://github.com/nf-core/tools/releases/tag/3.5.0).

As always, if you have any problems or run into any bugs, reach out on the [#tools slack channel](https://nfcore.slack.com/archives/CE5LG7WMB).

## Resolving conflicts on pipeline sync PRs

With this release, there are few merge conflicts, but they are still worth mentioning:

### Changing `Channel` to `channel`

As part of the new Nextflow strict syntax, we have changed all the mentions of `Channel` to `channel` (note the lower case `c`).
This might have introduced some merge conflicts in your `*.nf` files if you have modified them from the template,
for example in `subworkflows/local/utils_nfcore_$PIPELINE_pipeline/main.nf`.

- Resolution \*

Keep your changes, and afterwards change all mentions of `Channel` to `channel`.
