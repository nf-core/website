---
title: Resource optimization for all nf-core test profiles
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C06MD9L16UF
image: "/assets/images/events/2025/hackathon-barcelona/optimize_resources_config_profile_image.jpg"
image_alt: "Shouting cartoon figure with the caption 'Optimize all the things'"
leaders:
  flowuenne:
    name: Florian Wuennemann
    slack: "https://nfcore.slack.com/archives/DTZKT23D1"
---

## Project Aim

The `test` and `test_full` profiles in nf-core provide excellent tooling for testing pipelines before production deployment. However, these profiles currently over-allocate resources by default, leading to wasted compute and unnecessary costs during CI/CD runs and performance testing.

This project builds on previous optimization work documented at: https://github.com/FloWuenne/megatest-resource-optimization

## Goals

Since this project aims at optimizing all nf-core pipeline test profiles, it is probably more suited for experienced Nextflow / nf-core contributors.

1. Develop automated tooling to generate optimized resource configurations for each nf-core pipeline (Difficulty: Medium)
2. Integrate this automation with existing CI/CD workflows and AWS Megatests (Difficulty: Medium/Hard)
3. Create a system that automatically applies / adds optimized configurations after successful AWS Megatest runs (Difficulty: Hard)

The end result will be more efficient test profiles that reduce computational overhead while maintaining testing reliability.
