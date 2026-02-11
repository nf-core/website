---
title: nf-core/modules topics migration
category: components
slack: https://nfcore.slack.com/archives/C09LJTQQ3EY
intro_video: "https://www.youtube.com/watch?v=x_rJu-2M0Yo"
location: online
leaders:
  atrigila:
    name: Anabella Trigila
    slack: https://nfcore.slack.com/team/U03HEFYUW3H
---

This hackathon project focuses on updating
[nf-core/modules](https://github.com/nf-core/modules) to use modern Nextflow
features, specifically **topic channels**.

Topic channels simplify version collection across modules,
and are becoming the standard approach within nf-core.

Resources:

- [Nextflow documentation](https://nextflow.io/docs/latest/reference/channel.html#topic)
- [Adoption statistics dashboard](https://nf-core-stats.netlify.app/code/container_conversion/#version-topics-adoption-over-time)

---

## Goal

Migrate all nf-core modules to use topic channels.

This work is organized so each participant updates one module at a time,
making it suitable as a beginner-friendly contribution.

---

## What participants will do

Each contributor will:

1. Pick one module that has not yet been migrated.
2. Update it to use topic channels following official guidelines.
3. Run tests and lint checks locally.
4. Open a Pull Request for review.

Each migration corresponds to one issue in the modules repository.

---

## Tasks

### Migrate a module to topic channels

1. Choose an available module from the [migration tracker](https://github.com/nf-core/modules/issues/9978)

2. Assign the issue to yourself to avoid duplicate work.

3. Follow the [official migration guide](https://nf-co.re/docs/tutorials/migrate_to_topics/update_modules)

4. Update the module by:
   - Replacing legacy version collection logic
   - Using topic channels for version reporting
   - Updating tests if necessary

5. Run module tests and lint checks locally.

6. Open a Pull Request referencing the migration issue.

7. Address CI or review feedback until merged.

After finishing, contributors are encouraged to migrate additional modules.

---

## Recommended preparation

Participants should ideally have:

- Basic familiarity with Git and GitHub (forking repositories, creating branches, and opening pull requests)
- Basic knowledge of Nextflow
- Familiarity with nf-core modules

The following training material is recommended:

- [Hello nextflow](https://training.nextflow.io/latest/hello_nextflow/)
- [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/)

---
