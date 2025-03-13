---
title: Processes are retrying
shortTitle: Retries
---

## Processes are retrying

### Why did processes report an error but then retry?

One of the nice things about Nextflow is it offers the ability to retry processes if they encounter an error and fail with certain bash [exit status](https://en.wikipedia.org/wiki/Exit_status) or codes.
Some of these errors have common causes, which allows us to provide solutions to these problems.

A common issue is a tool requiring more memory than is (initially) made available to it based on the default memory specifications set in the pipeline.
Such errors (out of memory, or OOM) are often identified by exit code `104`, and those falling between `130` - `145`.

Therefore all nf-core pipelines [by default will retry](https://github.com/nf-core/tools/blob/930ece572bf23b68c7a7c5259e918a878ba6499e/nf_core/pipeline-template/conf/base.config#L18) a process if it hits one of those exit codes, but requesting more resources (memory, CPUs, and time) for the re-submitted job.

All other exit codes will cause the pipeline to fail immediately, and will not be retried.

However, some pipelines may extend this list or provide different retry conditions based on the behaviour of the specific tools used in the pipeline.
