---
title: Commit strategy
subtitle: Best practices for managing commits and merges
weight: 5
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Merge commit preference

We prefer to use merge commits when merging pull requests.
This approach helps avoid merge conflicts when multiple people are working with overlapping feature branches.

:::info{title="Rationale" collapse}
While merge commits create a more verbose git history, they make it much easier to manage parallel development.
When multiple contributors are working on related features or pulling from overlapping branches, merge commits preserve the individual branch histories and reduce the likelihood of conflicts.
:::

## Commit history philosophy

We prefer verbose commit histories but easy merges.

This means:

- Individual commits within a pull request can be small and frequent
- The full development history is preserved
- Each merge creates a clear record of when features were integrated
- Conflicts are minimized through the merge commit strategy

## Feature branch cleanup

Feature branches SHOULD be immediately deleted after merge.

The easiest way to ensure this happens is to enable the automatic deletion feature in your GitHub repository settings.
This keeps the repository clean and makes it clear which work is complete.

:::note
Squashing commits in a PR before merging is acceptable if preferred by the individuals working on that pull request.
The key is to choose an approach that works for your team while maintaining a clear history.
:::
