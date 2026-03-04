---
title: Subworkflow parameters
subtitle: Use subworkflow parameters
weight: 4
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Usage of parameters

Subworkflows MUST NOT assume that named `params` defined in the parent workflow will be passed to the subworkflow.
This allows developers to name their parameters as needed.
Use additional `input` value channels for such scenarios.
