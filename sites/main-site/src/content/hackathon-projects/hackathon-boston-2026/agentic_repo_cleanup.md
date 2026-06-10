---
title: Agentic Swarm Repo Cleanup
category: components
slack: https://nfcore.slack.com/archives/C0B1AJFK9MW
image: /assets/images/events/2026/hackathon-boston/agentic_repo_cleanup.png
image_alt: An AI swarm of agents collaboratively cleaning up a repository
leaders:
  HReed1:
    name: Harrison Reed
    slack: https://nfcore.slack.com/team/U0AQS30DMQT
---

Scaling tedious-but-necessary repository maintenance using AI agentic systems.

**nf-core/modules** contains hundreds of modules, many of which carry accumulated
technical debt: missing stub blocks, inconsistent test coverage, outdated
conventions, and absent documentation. These are the tasks that human developers
rarely _want_ to do — but they're essential for the health of the ecosystem.

This project demonstrates how **agentic AI swarms** can systematically tackle
this maintenance at scale, while keeping a human firmly in the loop for review
and submission.

### Track Record

Before the hackathon, this approach already delivered results:

- **48 modules** identified as missing stub blocks ([#4570](https://github.com/nf-core/modules/issues/4570))
- **40 modules** stubbed and validated using automated AST analysis and nf-test snapshot verification
- **10 atomic PRs** submitted ([#11349](https://github.com/nf-core/modules/pull/11349)–[#11358](https://github.com/nf-core/modules/pull/11358)), split by domain category for focused review
- Every stub follows the canonical `echo | gzip >` convention ([#5409](https://github.com/nf-core/modules/pull/5409))

### Resources

The methodology, tooling, and lessons learned are documented in detail:

- 🌐 [hvrinformatics.com](https://hvrinformatics.com) — Engineering blog covering agentic bioinformatics development
- 📝 [From Sandbox to Open Source: Deploying an Autonomous Swarm Against nf-core/modules](https://hvrinformatics.com/blog/nfcore-hackathon-stubs) — Deep dive into the stub generation pipeline
- 📝 [From 1 PR to 10: How an AI Agent Split 40 Modules Into Atomic Pull Requests](https://hvrinformatics.com/blog/nfcore-atomic-pr-strategy) — How the agent decomposed a 40-module PR into reviewable units
- 🛠️ [hvr-agentic-os](https://github.com/HReed1/hvr-agentic-os) — Open-source agentic swarm framework powering this project

---

## Goal

Use agentic AI systems to identify and resolve systematic maintenance gaps
across nf-core/modules — starting from the proven stub block foundation and
expanding into additional cleanup workstreams during the hackathon.

---

## Workstreams

### 1. Stub Block Completion (Proven)

Continue the [#4570](https://github.com/nf-core/modules/issues/4570) initiative.
Identify any remaining modules missing stubs, generate convention-compliant stub
blocks, validate with nf-test, and submit atomic PRs.

### 2. Additional Cleanup Targets (Hackathon Scope)

Potential areas to tackle collaboratively during the event:

- **nf-test coverage gaps** — Modules with incomplete or missing test snapshots
- **Topic channel migration** — Updating modules to use modern Nextflow topic channels
- **Linting compliance** — Resolving persistent `nf-core modules lint` warnings at scale
- **Documentation standardization** — Ensuring consistent `meta.yml` and module documentation

The exact scope for workstreams 2+ will be refined at the hackathon based on
community priorities and participant interest.

---

## What Participants Will Do

1. **Pick a workstream** — Choose from stub completion, test coverage, topic migration, or linting cleanup.
2. **Use the agentic workflow** — or contribute manually, both are welcome.
3. **Validate locally** — Run `nf-core modules lint` and `nf-test` before submission.
4. **Submit atomic PRs** — One module or small batch per PR for focused review.

This project is suitable for **all experience levels**. Each task is self-contained
and follows established nf-core contribution patterns.

---

## Recommended Preparation

- Basic familiarity with Git and GitHub
- Basic knowledge of Nextflow and nf-core module structure
- Review the [nf-core modules contribution guide](https://nf-co.re/docs/contributing/modules)
- Optional: Read the [stub block conventions from #5409](https://github.com/nf-core/modules/pull/5409)
- Optional: Explore [hvr-agentic-os](https://github.com/HReed1/hvr-agentic-os) to understand the agentic workflow
