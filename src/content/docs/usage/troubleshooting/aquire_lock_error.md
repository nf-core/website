---
title: Unable to acquire lock
subtitle: How to troubleshoot common mistakes and issues
weight: 70
---

# Unable to acquire lock error

Errors like the following:

```bash
Unable to acquire lock on session with ID 84333844-66e3-4846-a664-b446d070f775
```

Normally suggest a previous Nextflow run (on the same folder) was not cleanly killed by a user (e.g. using ctrl + z to hard kill a crashed run).

To fix this, you must clean the entirety of the run's `work/` directory e.g. with `rm -r work/` and re-running from scratch.

`ctrl +z` is **not** a recommended way of killing a Nextflow job. Runs that take a long time to fail are often still running because other job submissions are still running. Nextflow will normally wait for those processes to complete before cleaning shutting down the run (to allow rerunning of a run with `-resume`). `ctrl + c` is much safer as it will tell Nextflow to stop earlier but cleanly.
