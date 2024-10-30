---
title: PR review guidelines
subtitle: Best practices to review contributions.
---

## Get reviewers and feedback

- Authors/Reviewers can request help for review at nf-core Slack `#request-review` channel.
- We encourage collaboration in PR, other contributors are welcome to fix linting errors (you can mention relevant bot commands) or make improvements to the code if necessary.

## Merging conditions

Authors/Reviewers should merge after one positive review and with no obvious open questions.

### Exceptions

- Exception: merges to master for releases requires **two** reviews.
- Exception: pseudo PR reviews for a first release needs **at least one review from core or maintainer team**.

### Major changes

- If there is a major change in any type of PR (from modules to pipelines) then it's ok to request a second review - either the author or reviewer can do this.
- If the reviewer feels uneasy about a remaining open question then they shouldn't leave only a comment or approval, but rather a request-changes.

## Scope of the review process

- **request changes** reviews can be [dismissed](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/dismissing-a-pull-request-review) by the PR author if considered out of date and resolved.
- PRs that have approving reviews but remaining open questions can be merged by the author after addressing the questions to a common-sense level of satisfaction.
- If a PR has comments from one reviewer with open questions, and requests to re-review but does not perform the re-review, the PR may be merged if the first reviewer does not finish the review after 3 months if the PR gets an independent approval from someone else.

## Commit strategy

We prefer to use merge commits to merge. This is to avoid merge conflicts when multiple people are pulling with overlapping feature branches.

- We prefer verbose commit histories but easy merges.
- Feature branches should be immediately deleted after merge. Ideally by selecting the new feature on GitHub repo settings
- Note: squashing commits in a PR and merging after is fine, that is up to the individuals.
