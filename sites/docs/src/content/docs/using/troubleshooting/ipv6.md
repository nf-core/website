---
title: Network unreachable on IPv6 machines
subtitle: How to troubleshoot common mistakes and issues
weight: 90
---

### The pipeline is not getting pulled eventhough internet connection is stable

When running Nextflow pipelines on a IPv6 machine, you might encounter network errors even with a working internet connection. Below is a guide to troubleshoot and resolve these issues, particularly when using Docker and IPv6.

You might encounter the following error when running a Nextflow pipeline on a IPv6 machine:

```console
nextflow run nf-core/fetchngs -r 1.12.0 -profile test,docker --outdir test

WARN: Cannot read project manifest -- Cause: Network is unreachable (connect failed)
Cannot find `nf-core/fetchngs` -- Make sure exists a GitHub repository at this address `https://github.com/nf-core/fetchngs`
```

Your machine uses IPv6, but Nextflow or Java defaults to IPv4, causing network issues. Set the Java network stack to prefer IPv6 by exporting this variable before running your pipeline or globally add it to your `~/.bashrc`:

```bash
export NXF_OPTS="-Djava.net.preferIPv4Stack=false -Djava.net.preferIPv6Addresses=true"
```

### A process crashes due to network issues within the Docker container

If a model tries to fetch information of download data by using a Docker contaner, this can also lead to network issues on a IPv6 machine, since containers are isolated and you usually don't have network connectivity inside of them. You can configure Docker to use IPv6 for the default bridge network, not IPv4. See the official documentation [here](https://docs.docker.com/config/daemon/ipv6/#use-ipv6-for-the-default-bridge-network).

If this did not resolve your issue, try the `--network host` flag when running the container:

Create a `custom.config` file that you want to add to your nextflow run via `-c`:

```groovy title="custom.config
process {
    withName: YOUR_MODULE {
        containerOptions = <boilerplate docker.runOptions/containerOptions> + '--network host'
```

Then try running your pipeline again:

```bash
nextflow run nf-core/fetchngs -r 1.12.0 -profile test,docker --outdir test -c custom.config
```
