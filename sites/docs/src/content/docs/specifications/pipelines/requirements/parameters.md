---
title: Standardised parameters
subtitle: Strive to have standardised usage.
menu:
  main:
    weight: 100
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines SHOULD use the same command line option names as other pipelines for comparable options where possible.
For example, `--input` and `--genome`.

Parameter names SHOULD work in a similar way across pipelines.
For example, `--input` typically takes a `.csv` sample sheet file where appropriate.

> We're planning to build a tool that lists every parameter used by every pipeline, so that you can check for existing parameters with similar names
> You can track the progress of this feature request here: [nf-core/website#1251](https://github.com/nf-core/website/issues/1251)
