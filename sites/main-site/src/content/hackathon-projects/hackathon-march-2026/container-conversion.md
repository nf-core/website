---
title: Seqera Containers migration
category: components
slack: https://nfcore.slack.com/archives/C09V4RAHV0X
location: online
intro_video: ""
leaders:
  mashehu:
    name: Matthias Hörtenhuber
    slack: https://nfcore.slack.com/team/UQZG2UCF3
---

This project continues the ongoing migration of nf-core modules from
[BioContainers](https://biocontainers.pro/) to
[Seqera Containers](https://seqera.io/containers/), built using
[Wave](https://seqera.io/wave/).

The migration enables automatic multi-architecture container builds
(`linux/amd64` and `linux/arm64`), Conda lock files for improved
reproducibility, and prepares the modules for fully automated update
pipeline via Renovate and GitHub Actions.

Resources:

- [Blog: Migration from Biocontainers to Seqera Containers — Part 2](https://nf-co.re/blog/2024/seqera-containers-part-2)
- [nf-core/tools PR #3954](https://github.com/nf-core/tools/pull/3954)
- [nf-core/tools PR #3955](https://github.com/nf-core/tools/pull/3955)
- [Bulk module update tracker (nf-core/modules #6698)](https://github.com/nf-core/modules/issues/6698)
- [Container conversion stats dashboard](https://nf-core-stats.netlify.app/code/container_conversion/)

---

## Goal

Migrate nf-core modules and pipelines to use Seqera Containers:

- **Modules**: build new containers with `nf-core modules containers create`, and let the tool update `meta.yml` and `main.nf` automatically.
- **Pipelines**: update installed modules and let `nf-core modules update` auto-generate the per-platform container config files in `conf/`.
- **tools**: Add missing functionality for the above goals and
  fix any related bugs.

---

## Tasks

### Convert a module to Seqera Containers

1. Claim an open module by opening a PR and assigning yourself to it.

2. Run the following `nf-core` command to build containers via the Wave API and
   update `meta.yml` and `main.nf` automatically:

   ```bash
   nf-core modules containers create <module_name>
   ```

   This builds Docker and Singularity images for both `linux/amd64` and
   `linux/arm64`, writes the image URLs and build IDs into `meta.yml`, and
   updates the `container` string in `main.nf`.

   :::info
   Depending on the current build work load and size of the container,
   this can take several minutes.
   :::

   :::note
   Set `TOWER_ACCESS_TOKEN` to your Seqera Platform token to avoid
   strict Wave API rate limits.
   :::

3. Verify the changes look correct, run tests, and open a Pull Request.

---

### Update a pipeline to use Seqera Containers

Once modules have been converted, pipelines can adopt the new containers by
updating their installed or local modules. The new tooling then auto-generates
per-platform container config files under `conf/` using `nextflow inspect`.

1. Inside a pipeline repository, update one or more modules:

   ```bash
   nf-core modules update <module_name>
   # or update all modules at once:
   nf-core modules update
   ```

2. The tool automatically runs `nextflow inspect` and writes config files for
   all supported platforms:

   :::note
   This requires Nextflow >= 25.04.4.
   :::

3. Review the generated config files, run the pipeline tests, and open a Pull
   Request.

---

## Recommended preparation

Participants should ideally have:

- Basic familiarity with Git and GitHub (forking, branching, pull requests)
- Basic knowledge of Nextflow and nf-core modules or python
