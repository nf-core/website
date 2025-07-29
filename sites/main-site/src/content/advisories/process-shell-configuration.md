---
title: .command.run Permission denied
subtitle: Nextflow versions after 24.10.6 are incompatible with certain nf-core pipelines published in 2024 and early 2025.
category:
  - pipelines
  - configuration
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
  - name: ampliseq
    versions:
      - 2.12.0
  - name: airrflow
    versions:
      - 4.2.0
  - name: bacass
    versions:
      - 2.4.0
  - name: crisprseq
    versions:
      - 2.3.0
  - name: demo
    versions:
      - 1.0.1
  - name: demultiplex
    versions:
      - 1.5.3
      - 1.5.4
  - name: denovotranscript
    versions:
      - 1.1.0
  - name: detaxizer
    versions:
      - 1.1.0
  - name: fastquorum
    versions:
      - 1.1.0
  - name: mag
    versions:
      - 3.2.0
      - 3.2.1
      - 3.3.0
  - name: metapep
    versions:
      - 1.0.0
  - name: methylseq
    versions:
      - 2.7.0
      - 2.7.1
      - 3.0.0
  - name: nanostring
    versions:
      - 1.3.1
  - name: rnaseq
    versions:
      - 3.16.1
      - 3.17.0
      - 3.18.0
  - name: pairgenomealign
    versions:
      - 1.1.1
  - name: phaseimpute
    versions:
      - 1.0.0
  - name: rangeland
    versions:
      - 1.0.0
  - name: sarek
    versions:
      - 3.5.0
  - name: scrnaseq
    versions:
      - 3.0.0
  - name: smrnaseq
    versions:
      - 2.4.0
  - name: taxprofiler
    versions:
      - 1.2.1
      - 1.2.2
modules:
subworkflows:
configuration:
nextflowVersions:
  - ">=25.04.0 <25.04.3"
nextflowExecutors:
softwareDependencies:
references:
  - title: Correction in nf-core pipeline template
    description: Pull-request that fixed the underlying issue in the nf-core pipeline template
    url: https://github.com/nf-core/tools/pull/3416
  - title: Announcement regarding the issue on nf-core Slack
    description: Detailed report by Phil Ewels regarding the matter and suggestions for resolution
    url: https://nfcore.slack.com/archives/CE6P95170/p1747302440656339
---

# Issue

After Nextflow 25.04 was released, many users reported getting `.command.run: Permission denied errors` immediately after launching recent pipeline versions of nf-core pipelines.
This is an obscure error that emerged due to an unfortunately invalid configuration syntax in the nf-core pipeline template, that no longer works due to a stricter configuration parsing in the most recent Nextflow versions.

# Resolution

If you're seeing this above error you have a few options. In order of difficulty:

- Use an older version of Nextflow (24.10.6 or earlier)
- Use a Nextflow config file that overwrites the invalid config syntax locally. Save this in a config file

  ```bash
    // Set bash options
  process.shell = [
      "bash",
      "-C",         // No clobber - prevent output redirection from overwriting files.
      "-e",         // Exit if a tool returns a non-zero status/exit code
      "-u",         // Treat unset variables and parameters as an error
      "-o",         // Returns the status of the last command to exit..
      "pipefail"    //   ..with a non-zero status or zero if all successfully execute
  ]
  ```

  and use with `-c`, eg. `nextflow run nf-core/rnaseq -c nf_config_fix.config`.

- Use an older version of the pipeline that has nf-core/tools template < 3.0.0, thus was published before 2024-10-08.
  The issue was fixed and published in the template 3.1.2 on 2025-01-27, so pipeline versions released after that day are mostly not affected either, unless they were published with an older template.
  You can fine the template version that a pipeline uses in a hidden file called `.nf-core.yml` in the pipeline repository.
  Additionally, some pipelines already incorporate a manual fix, so even though released during that period, pixelator v1.4.0 and references 0.1 already incorporate the required change.
