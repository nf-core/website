---
title: Advisories
subtitle: Learn how to create and publish advisories for nf-core components
weight: 3
---

Advisories serve as structured, long-lived notices for significant technical issues affecting nf-core components.
They prioritise clarity and practical value, helping users avoid common pitfalls.

<!-- TODO: Add image of advisory somewhere on this page -->

## What are advisories?

Advisories address problems substantial enough to consume hours of user troubleshooting. They provide structured information about:

- Bugs causing incorrect results in modules or pipelines
- Known regressions with specific parameter combinations
- Version incompatibilities between pipelines and Nextflow
- Security vulnerabilities in dependencies or containers
- Executor-specific problems requiring special handling

:::note
Advisories differ from regular issues or documentation updates.
They focus on significant problems that require structured, searchable, long-term visibility for the community.
:::

## Who can publish advisories?

Anyone can create an advisory.
Whilst maintainers typically author advisories, the process welcomes contributions from any community member who identifies a worthy issue.

## Publishing an advisory

To publish and advisory:

1. Fork and add content:
   1. Fork the [nf-core website repository](https://github.com/nf-core/website)
   1. Add a new Markdown or MDX file to `sites/main-site/src/content/advisory`
1. Add frontmatter:
   - Reference existing advisories as templates for structure and formatting
   - Use the following as a template:

     ```yaml
     ---
     # Required fields
     title: "Advisory Title"
     subtitle: "Brief description of the issue"
     category:
       "pipelines" # Which part of the nf-core ecosystem this advisory affects
       # Options: ["pipelines", "modules", "subworkflows", "configuration"]
       # Can be single value or array for multi-category issues
     type:
       "known_regression" # What kind of issue this is - helps users understand impact
       # Options: ["known_regression", "incompatibility", "security",
       #          "performance", "data_corruption", "other"]
       # Can be single value or array for issues with multiple aspects
     severity:
       "high" # How serious this issue is for users
       # Options: ["low", "medium", "high", "critical"]
       # Note: "critical" is only allowed for security issues
     publishedDate: "2024-01-15" # When this advisory was published (YYYY-MM-DD format)

     # Optional reporter information - who discovered and reported this issue
     reporter: # Can be null, array of usernames, or array of objects with details
       - "username" # Simple GitHub username
       - name: "Full Name" # Object with full name and GitHub username
         github: "username"

     # Optional reviewer information - who reviewed and validated this advisory
     reviewer: # Array of usernames or objects with reviewer details
       - "reviewer-username"

     # Category-specific fields - REQUIRED if the corresponding category is specified above
     pipelines: # List of affected pipelines (required if category includes "pipelines")
       - "pipeline-name" # Simple list of pipeline names
       # OR specify affected versions:
       - name: "pipeline-name" # Pipeline name with specific version information
         versions: ["1.0.0", "1.1.0"] # Semantic version numbers of affected releases

     modules: # List of affected modules (required if category includes "modules")
       - "module_name"
       - "another/module"

     subworkflows: # List of affected subworkflows (required if category includes "subworkflows")
       - "subworkflow/name"

     configuration: # Configuration aspects affected (required if category includes "configuration")
       - "config-item"

     # Optional details
     nextflowVersions: # Specific Nextflow versions that exhibit this issue
       - "23.04.0" # Use semantic versioning format
       - "23.10.1"

     nextflowExecutors: # Workflow execution environments where this issue occurs
       - "SLURM" # Specific executor names
       - "AWS Batch"
       - "Local"

     softwareDependencies: # Container systems or package managers affected by this issue
       - "Docker" # Simple list of affected systems
       # OR specify affected versions:
       - name: "Singularity" # System name with version details
         versions: ["3.8.0", "3.9.0"] # Specific versions that have the issue
       # Available: ["Apptainer", "Charliecloud", "Docker", "Podman", "Sarus",
       #            "Shifter", "Singularity", "Conda", "Spack", "Wave"]

     # Optional references - links to related information, bug reports, documentation
     references:
       - title: "GitHub Issue" # Short title describing what this link is
         description: "Original bug report" # Longer explanation of what you'll find at this link
         url: "https://github.com/nf-core/pipeline/issues/123" # The actual URL
       - title: "Slack discussion"
         description: "Original reporting on the nf-core slack"
         url: "https://nfcore.slack.com/archives/C03EZ806PFT/p1730391850337429"
     ---
     ```

1. Submit your advisory
   1. Create a pull request to the `nf-core/website` repository
   1. Request review from maintainers
   1. Address any feedback and update your advisory as needed
