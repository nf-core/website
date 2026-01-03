---
title: Use nf-core git branches
subtitle: Use `master|main`, `dev` and `TEMPLATE`.
menu:
  main:
    weight: 160
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

The latest stable release MUST be on the main `main` or `master` branch.
No additional changes SHOULD be pushed to master after each release.

The main development code MUST be kept in a branch called `dev`.
The nf-core `dev` branch SHOULD be set as the default branch up until the first release.

For minor bugfixes a `patch` branch MAY be used and merged directly into `master|main`, leaving `dev` for continued development work.

The `TEMPLATE` branch MUST only contain vanilla nf-core template code.
It is used for automated synchronisation of template updates.
