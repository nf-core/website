---
title: Bundled documentation
subtitle: Pipeline documentation must be hosted on the nf-core website
menu:
  main:
    weight: 60
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

All documentation MUST be bundled with the pipeline code in the main repository, within a directory called `docs`.

Documentation MUST be hosted only on the nf-core website and GitHub.
Hosting the documentation at a second location (such as custom readthedocs website, or GitHub pages) MUST NOT be permitted.

:::info{title="Rationale" collapse}
This ensures that users of nf-core pipelines can always intuitively find the documentation for all nf-core pipelines in the same way, providing a consistent user experience.
:::

Documentation MUST include at least the following files:

- `README.md`
- `docs/usage.md`
- `docs/output.md`

Additional pages (e.g., tutorials, FAQs) MAY be added and will be automatically rendered on the nf-core website pipeline page.
