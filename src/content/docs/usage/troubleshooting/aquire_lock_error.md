---
title: Unable to acquire lock
subtitle: How to troubleshoot common mistakes and issues
shortTitle: Unable to acquire lock
weight: 70
---

# Unable to acquire lock error

Errors like the following:

```bash
Unable to acquire lock on session with ID 84333844-66e3-4846-a664-b446d070f775
```

This error suggests that a user did not clealy kill a previous Nextflow run in the same launch directory.

To fix this, delete the [workDir](https://www.nextflow.io/docs/latest/config.html#miscellaneous), which is called `work` by default.

:::warning
`ctrl +z` is **not** a recommended way of killing a Nextflow job. Runs that take a long time to fail are often still running because other job submissions are still running. Nextflow will normally wait for those processes to complete before cleaning shutting down the run (to allow rerunning of a run with `-resume`). `ctrl + c` is much safer as it will tell Nextflow to stop earlier but cleanly.
:::
