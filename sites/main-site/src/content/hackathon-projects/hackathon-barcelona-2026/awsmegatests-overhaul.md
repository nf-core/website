---
title: AWS megatests overhaul
category: tooling
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/megatests.jpg"
image_alt: "Marie Kondo loves AWS megatests"
slack: https://nfcore.slack.com/archives/CTY8L26C8
leaders:
  FriederikeHanssen:
    name: Friederike Hanssen
    slack: https://nfcore.slack.com/team/UPC8CHDKQ
---

## Goal

Overhaul the nf-core AWS megatests: the full-size test runs that every pipeline executes on release.

## Tasks

_Tasks will be added closer to the event. Check back for updates._

- Update CE definitions to be more stable for on-demand fallbacks
- Port CE definitions to Terraform
- Change user management to Teams
- Add 'Write' permission to launch users, so they can add pipelines to the launchpad. Consider adding pipelines to the launchpad automatically and manage new releases with pipeline versions
- Update platform action to use the latest release:
  - Make use of new features, such as wait on submission
- Write and update user docs + docs for the core team:
  - How to switch CEs on failures
  - Most common work arounds
  - What to do when the action fails
- Wipe CEs on AWS account and remove all invalid CEs
- Consider if we can get a dashboard for which pipeline has full size data available on which release

:::note
Keep an eye on this page for more updates.
:::
