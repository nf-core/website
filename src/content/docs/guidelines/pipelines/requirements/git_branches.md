---
title: Use nf-core git branches
subtitle: Use `master`, `dev` and `TEMPLATE`.
menu:
  main:
    weight: 160
---

The latest stable release should be on the main `master` branch.
No additional changes should be pushed to master after each release.

The main development code should be kept in a branch called `dev`.
The nf-core `dev` branch should be set as the default branch up until the first release.

For minor bugfixes a `patch` branch may be used and merged directly into `master`, leaving `dev` for continued development work.

The `TEMPLATE` branch should only contain vanilla nf-core template code.
It is used for automated synchronisation of template updates.
