---
title: Update pipelines to use topic channels
description: How to update pipelines to use the new topic channels feature
---

Adapting pipelines to use topic channels for version outputs is a straightforward process explained in the steps below:

1. Update the template of the pipeline to the latest version that includes support for topic channels. This is supported starting from version 3.5.0.

2. Pull the latest changes made to modules using `nf-core modules update`.

3. Look for issues in the pipeline code where `versions` is an invalid output of a process. Remove these outputs since these are now handled by the topic channels.
