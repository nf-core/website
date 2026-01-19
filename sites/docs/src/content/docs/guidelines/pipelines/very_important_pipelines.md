# Very Important Pipelines

## Overview

Widely adopted nf-core pipelines play a critical role in the ecosystem. They are widely used by the community and often act as reference for the development of the rest of the pipelines.

This document defines a tier for such pipelines, called **Very Important Pipelines**, with elevated expectations around maintenance, governance, testing, and documentation. The intent is to formalize best practices that many popular pipelines already follow and to provide a clear model for long term, sustainable pipeline development.

## Scope

Very Important Pipelines are standard nf-core pipelines that have reached high community adoption and visibility. All existing nf-core requirements continue to apply. This tier adds additional recommendations appropriate for pipelines with a large user and contributor base.

This tier is additive, not exclusive.

## Identification of Very Important Pipelines

There are no strict thresholds. Selection is based on an overall assessment by the nf-core core team, considering factors such as:

* Number of pull requests and contributors
* Repository stars and forks
* Pipeline age and stability
* Issue volume and responsiveness
* Website traffic and usage indicators

Selection is expected to evolve organically as pipelines grow and mature.

## Team and Governance Recommendations

Very Important Pipelines should demonstrate shared ownership and active governance.

Recommended practices:

* A team of active co-developers
* A dedicated GitHub team
* Defined code owners using CODEOWNERS
* Regular development meetings, open to the community
* Use of milestones to plan and track releases
* A developer focused Slack channel, separate from user support

## Development and Maintenance Recommendations

Very Important Pipelines should maintain a high standard of technical quality and currency.

Recommended practices:

* Regular review of default configuration usage
* Pipeline template and tool versions kept within two releases of current
* Automated testing coverage above eighty percent
* Clear issue advisories when known problems are identified
* Participation in an nf-core bytesize talk

## Documentation and User Experience

Very Important Pipelines should set the standard for clarity and usability.

Recommended practices:

* A pipeline logo and associated visual assets
* A pipeline metro map
* Extended documentation beyond the minimum requirements
* All outputs fully described in the documentation
* Practical usage tutorials for common workflows

## Implementation

The nf-core core team will:

* Integrate this guidance into the nf-core website
* Identify an initial set of Very Important Pipelines
* Reach out to pipeline teams to support adoption of these recommendations

Designation is expected to be collaborative rather than prescriptive.

## Branding

Very Important Pipelines will display an additional badge in the repository README and on the nf-core website to indicate their status.

## Goals

* Increase visibility of widely adopted nf-core pipelines
* Promote sustainable development practices for high-impact pipelines
* Provide a clear model for emerging pipelines to follow

## Non-Goals

* This tier does not devalue or discourage smaller or newer pipelines
* It does not replace existing nf-core standards
* Selection criteria are indicative and may change over time
