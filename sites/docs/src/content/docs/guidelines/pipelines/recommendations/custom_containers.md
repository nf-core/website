---
title: Custom containers
subtitle: What to do with custom containers that are hosted on docker.io or ghcr.io
menu:
  main:
    weight: 200
---

If a pipeline cannot be constructed with bioconda and therefore cannot be found on biocontainers, please contact `@core-team` on slack to discuss the best way to proceed.
The current solution to ensure future reproducibitliy is mirror the container and host it on `quay.io` in the nf-core organisation.
For that the following steps are necessary, and should be done by a core-team member.

```bash
docker pull <ORGANISATION>/<TOOL>:<TAG>
docker tag <ORGANISATION>/<TOOL>:<TAG> quay.io/nf-core/TOOL:<TAG>
docker push quay.io/nf-core/TOOL:<TAG>
```

Once this is done, please go on to the quay.io website, and check the newly updated repository's settings (that's the quay.io name for the container) and make sure that its visibility is public.
