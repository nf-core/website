---
title: PR review guidelines
subtitle: Best practices to review contributions.
markdownPlugin: addNumbersToHeadings
---

The key words “MUST”, “MUST NOT”, “SHOULD”, etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Actively request reviewers and feedback

- Authors/Reviewers SHOULD request help for review at nf-core Slack `#request-review` channel for general reviews (please use [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) for pipeline releases)
- We encourage collaboration in PR, other contributors are welcome to fix linting errors (you can mention relevant bot commands) or make improvements to the code if necessary.

## Pull Request Rapid Merge After Review

Authors/Reviewers SHOULD merge after one positive review and with no obvious open questions.

### Exceptions

- Pull requests to master for releases MUST receive **two** positive reviews before merging.
- (Pseudo) pull requests to master for [_first_ releases](/docs/tutorials/adding_a_pipeline/first_release) MUST receive **at least one review from core or maintainer team**.

### Addressing Major Change Reviews

- If there is a major change in any type of PR during review (from modules to pipelines) then developers SHOULD request a second review - either the author or reviewer can do this.
- If the reviewer has strong reservations about proposed changes, a reviewer SHOULD NOT leave only a comment or approval, but SHOULD leave a request-changes review.

## Review Process Scope

- **Request changes** reviews MAY be [dismissed](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/dismissing-a-pull-request-review) by the PR author if considered out of date and resolved.
- PRs that have approving reviews but remaining open questions MAY be merged by the author after addressing the questions to a common-sense level of satisfaction.
- If a PR has comments from one reviewer with open questions, and they request to re-review but does not perform the re-review, the PR MAY be merged if the first reviewer does not finish the review after 3 months if the PR gets an independent approval from someone else.

## Commit strategy

We prefer to use merge commits to merge. This is to avoid merge conflicts when multiple people are pulling with overlapping feature branches.

- We prefer verbose commit histories but easy merges.
- Feature branches SHOULD be immediately deleted after merge. Ideally by selecting the new feature on GitHub repo settings
- Note: squashing commits in a PR and merging after is fine, that is up to the individuals.
