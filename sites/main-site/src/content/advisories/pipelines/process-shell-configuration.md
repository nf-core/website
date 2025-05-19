---
title: .command.run Permission denied
subtitle: Nextflow versions after 24.10.6 are incompatible with certain nf-core pipelines published in 2024 and early 2025.
category: pipelines
type: incompatibility
severity: high
publishedDate: "2025-05-15"
reporter:
 - alexnater
 - MatthiasZepper
reviewer:
  - ewels
  - MatthiasZepper
pipelines:
- name: rnaseq
  versions:
    - 3.16.1
    - 3.17.0
    - 3.18.0
modules:
subworkflows:
configuration:
nextflowVersions:
  - 25.04.0
  - 25.04.1
  - 25.04.2
nextflowExecutors:
softwareDependencies:
references:
  -   title: Correction in nf-core pipeline template
      description: Pull-request that fixed the underlying issue in the nf-core pipeline template
      url: https://github.com/nf-core/tools/pull/3416
  -   title: Announcement regarding the issue on nf-core Slack
      description: Detailed report by Phil Ewels regarding the matter and suggestions for resolution
      url: https://nfcore.slack.com/archives/CE6P95170/p1747302440656339
---

# Issue

After Nextflow 25.04 was released, many users reported getting `.command.run: Permission denied errors` immediately after launching recent pipeline versions of nf-core pipelines. This is an obscure error that emerged due to an unfortunately invalid configuration syntax in the nf-core pipeline template, that no longer works due to a stricter configuration parsing in the most recent Nextflow versions.

# Resolution

If you're seeing this above error you have a few options. In order of difficulty:

- Use an older version of Nextflow (24.10.6 or earlier)
- Use a Nextflow config file that overwrites the invalid config syntax locally (see thread for details :thread:)
- Use an older version of the pipeline that has nf-core/tools template < 3.0.0, thus was published before 2024-10-08. The issue was fixed and published in the template on 2025-01-27, so pipeline versions released after that day are presumably not affected either.
