---
title: Custom Docker images
subtitle: What to do with custom containers that are hosted on docker.io or ghcr.io
menu:
  main:
    weight: 200
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

If a pipeline cannot be constructed with bioconda and therefore cannot be found on biocontainers, please contact `@core-team` on Slack to discuss the best way to proceed.

Where possible, we prefer to mirror the container and host it on `quay.io` in the nf-core organisation to ensure future reproducibility.

:::info{title="Rationale" collapse}
Mirroring ensures that a container that the pipeline relies on will not be deleted by a third party in the future.
:::

To request the mirroring, please request via the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel.

A core-team member will run the following commands:

```bash
docker pull <organisation>/<tool>:<tag>
docker tag <organisation>/<tool>:<tag> quay.io/nf-core/tool:<tag>
docker push quay.io/nf-core/tool:<tag>
```

Once this is done, please go on to the quay.io website, and check the newly updated repository's settings (that's the quay.io name for the container) and make sure that its visibility is public.
