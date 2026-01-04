---
title: Rapid merge after review
subtitle: Guidelines for merging pull requests after approval
weight: 3
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Pull request rapid merge after review

To maintain momentum and avoid unnecessary delays, authors and reviewers SHOULD merge pull requests promptly after receiving approval.

A pull request SHOULD be merged after:

- One positive review has been received, AND
- There are no obvious open questions or unresolved concerns

## Exceptions

### Release pull requests

Pull requests to the master branch for releases have stricter requirements:

- MUST receive **two** positive reviews before merging

### First releases

Pseudo pull requests to master for [first releases](/docs/tutorials/adding_a_pipeline/first_release) have special requirements:

- MUST receive **at least one review from core or maintainer team**

## Addressing major change reviews

When significant changes occur during the review process:

- If there is a major change in any type of PR (from modules to pipelines), developers SHOULD request a second review
- Either the author or reviewer can request this additional review

### Reviewer reservations

If a reviewer has strong reservations about proposed changes:

- The reviewer SHOULD NOT leave only a comment or approval
- The reviewer SHOULD leave a "request changes" review instead
