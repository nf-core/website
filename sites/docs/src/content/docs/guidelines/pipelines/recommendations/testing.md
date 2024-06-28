---
title: Testing
subtitle: Use nf-test to validate pipeline
weight: 210
---

Pipelines should use nf-test for valid pipeline testing.

All pipelines should support nf-test and use it for pipeline level tests. Additional tests for local subworkflows and modules is recommended. Modules and subworkflows from nf-core/modules should include nf-tests which can also be used within the pipeline.
