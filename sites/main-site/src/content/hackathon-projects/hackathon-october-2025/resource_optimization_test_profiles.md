---
title: Resource optimization for all nf-core test profiles
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C06MD9L16UF
image: "/assets/images/events/2025/hackathon-barcelona/optimize_resources_config_profile_image.jpg"
image_alt: "Optimize all the things"
leaders:
  flowuenne:
    name: Florian Wuennemann
    slack: "https://nfcore.slack.com/archives/DTZKT23D1"
---

## Project Aim

The `test` and `test_full` profiles in nf-core provide excellent tooling for testing pipelines before production deployment. However, these profiles currently over-allocate resources by default, leading to wasted compute and unnecessary costs during CI/CD runs and performance testing.

This project builds on previous optimization work documented at: https://github.com/FloWuenne/megatest-resource-optimization

## Goals

1. Develop automated tooling to generate optimized resource configurations for each nf-core pipeline
2. Integrate this automation with existing CI/CD workflows and AWS Megatests
3. Create a system that automatically opens pull requests with optimized configurations after successful test runs
4. Establish a standardized approach for resource optimization across all nf-core pipelines

The end result will be more efficient test profiles that reduce computational overhead while maintaining testing reliability.
