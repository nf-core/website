---
title: Maintenance and deprecation
subtitle: Follow config maintenance and deprecation guidelines
markdownPlugin: addNumbersToHeadings
weight: 1
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Maintenance

### Primary maintainer(s) responsibility

The names listed in the `config_profile_contact` SHOULD be responsible for maintaining the config as the primary maintainer.
The primary maintainer MAY be more than one person.

Primary maintainers SHOULD test any changes to configs on their infrastructure.
Primary maintainers SHOULD merge their own pull requests with changes.

If primary maintainers no longer have access to the infrastructure, they SHOULD make a best effort to transfer maintainership to another contact person who still has access.

### nf-core community

The nf-core community core- and maintainers teams MAY make minor modifications to configs on the nf-core/configs repository to ensure compatibility with Nextflow and the nf-core template.
The nf-core community core- and maintainers teams MUST NOT make major changes without permission and testing from the primary contact, for example, changing conditionals for selecting queues.
Should the primary maintainers become unresponsive, the nf-core community core- and maintainers teams MAY transfer maintainership to another, active, member of the institution willing to take the responsibility on.

## Deprecation

Deprecation of a config involves the deletion of all files related to that particular cluster (config, pipeline-specific configs, and documentation) from the nf-core/configs repository.

### Deprecation through change of primary maintainer

A primary maintainer MAY only deprecate a config if they are no longer able to maintain it and found no person to transfer maintainership to.

### Deprecation through decommissioning

Primary maintainers SHOULD deprecate a config if a particular infrastructure is decommissioned.

### Deprecation by the nf-core community

The nf-core community core- and maintainers teams MAY deprecate a config in certain circumstances.

The nf-core community core- and maintainers teams MAY only deprecate a config if problems with the config affect the entire repository of configs.

nf-core community core- and maintainers MUST follow the procedure below for deprecation:

1. Contact the lead maintainer about fixing and testing the config via GitHub, Slack, and email (if available).
2. Modify the config's `config_profile_description` with a large warning for users to contact the nf-core community.
3. Deprecate the config if no contact is made by the primary maintainer or users within two months.
