---
title: Custom Docker images
subtitle: What to do with custom containers that are hosted on docker.io or ghcr.io
menu:
  main:
    weight: 200
---

If a pipeline cannot be constructed with bioconda and therefore cannot be found on biocontainers, please contact `@core-team` on slack to discuss the best way to proceed.

If possible, we prefer to mirror the container and host it on `quay.io` in the nf-core organisation to ensure future reproducibility (ie. that a container that the pipeline relies on will not be deleted by a 3rd party in the future).

To request the mirroring, please request the mirroring via the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel

Then, a core-team member will run the following commands:

```bash
docker pull <ORGANISATION>/<TOOL>:<TAG>
docker tag <ORGANISATION>/<TOOL>:<TAG> quay.io/nf-core/TOOL:<TAG>
docker push quay.io/nf-core/TOOL:<TAG>
```

Once this is done, please go on to the quay.io website, and check the newly updated repository's settings (that's the quay.io name for the container) and make sure that its visibility is public.
