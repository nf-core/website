---
title: Agentic module development
category: components
slack: https://nfcore.slack.com/archives/CJRH30T6V
location: Barcelona
leaders:
  vagkaratzas:
    name: Evangelos Karatzas
    slack: https://nfcore.slack.com/team/U05LNHCFLCW
---

## Goal

Let agentic coding tools loose on the [nf-core/modules wishlist](https://github.com/nf-core/modules/issues?q=is%3Aissue%20state%3Aopen%20label%3Awishlist), with zero manual intervention, then have humans review and evaluate what comes out — so we learn where the agents actually fail and fix the instructions we give them. Then, of course, address any review points and merge to master.

## Description

Agentic coding tools (Claude Code, Codex, and anything else that supports plugins, skills and sub-agents) can already scaffold a module, write `nf-test`s and fill in a `meta.yml`.
What we do not know is _how well_, and _how much of that quality comes from the context we hand them_.

This project runs that as an experiment during the hackathon.
The group generates modules, and then trades and reviews them.

**Track A — generation.** Pick an open wishlist issue and hand it to an agent using exactly **one** of four setups:

| Method | Context given to the agent                                                                    |
| ------ | --------------------------------------------------------------------------------------------- |
| (i)    | The [`nf-core-module-dev`](https://github.com/vagkaratzas/nf-core-module-dev) plugin only        |
| (ii)   | The `AGENTS.md` of nf-core/modules only                                                         |
| (iii)  | Both the plugin and `AGENTS.md`                                                                 |
| (iv)   | Neither — a bare agent with only the issue text                                                 |

The agent is allowed to open the PR itself once it judges the issue is solved.
**No manual fixes, no "just one quick edit" before submitting** — a nudged PR is a contaminated data point.
The PR description must state which of the four methods was used, so reviewers and the final analysis can group results.

**Track B — review.** Other group members (using AI assistance or not, their call) review the incoming PRs against a shared rubric.

Scoring criteria, each 1–5 plus free-text notes:

- **PR readiness** — could this be merged, or does it need a round of fixes first?
- **Correctness** — does the module do what the tool actually does? Right arguments, right behaviour?
- **Fullness** — are all input files present and declared as they should be? Optional inputs, `meta` maps, all expected outputs?
- **Standards** — naming, `environment.yml`, `meta.yml` completeness, topic channels for versions, linting, snapshot quality
- **Tool difficulty** — how hard was this tool to wrap in the first place? A trivial CLI failing means something very different from a gnarly multi-output tool failing
- **Free-text notes** — the part we care most about: what specifically went wrong, and what context would have prevented it

## Tasks

Good project for anyone curious about agentic coding, and equally good for experienced module reviewers who have never touched an agent, since Track B can benefit from people who know what a correct module looks like.

:::info{title="Great first issues"}
This project is newcomer-friendly:
Track A just needs an agentic coding tool set up on your machine before you start.
Track B is reviewing a module PR against a rubric, and is one of the best ways to learn nf-core module standards.
:::

Specific tasks:

- Set up an agentic coding tool and, for Track A methods (i)/(iii), install the `nf-core-module-dev` plugin
- Claim a wishlist issue and a method on the shared tracking sheet, so we get a spread across all four methods and avoid duplicate work
- Run the agent hands-off and let it open the PR, with the method named in the PR body
- Record what you observed while it ran: retries, dead ends, where it got stuck
- Review incoming PRs against the rubric and post scores plus free-text notes
- Collate scores at the end and look for patterns per method and per criterion
- Turn the findings into concrete PRs against the plugin's skills/agents and against nf-core/modules `AGENTS.md`

Wishlist issues live here: [nf-core/modules wishlist](https://github.com/nf-core/modules/issues?q=is%3Aissue%20state%3Aopen%20label%3Awishlist).

If an agent fails in an interesting way that no criterion captures, write it up anyway — the failure modes we did not anticipate are the most useful output of this project.

## Making PRs

Module PRs go against the `master` branch of [nf-core/modules](https://github.com/nf-core/modules).
Plugin improvements go against [vagkaratzas/nf-core-module-dev](https://github.com/vagkaratzas/nf-core-module-dev), and `AGENTS.md` improvements against nf-core/modules.

Since this is an in-person project, grab Vangelis at the table, or post in the project Slack channel to get paired with a reviewer.
Note that agent-generated PRs will sit unmerged until they have been reviewed and scored, so please do not fast-track your own.
