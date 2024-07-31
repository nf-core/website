---
title: Bundled documentation
subtitle: Pipeline documentation must be hosted on the nf-core website
menu:
  main:
    weight: 60
---

All documentation must be bundled with the pipeline code in the main repository, within a directory called `docs`.

Documentation must _only_ be hosted on the nf-core website and GitHub.
Hosting the documentation at a second location (such as custom readthedocs website, or GitHub pages etc) is not allowed.
This is to ensure that users of nf-core pipelines can always intuitively find the documentation for all nf-core pipelines in the same way.

Documentation must include at least the following files:

- `README.md`
- `docs/usage.md`
- `docs/output.md`

Additional pages (e.g. tutorials, FAQs) can be added and will be automatically rendered on the nf-core website pipeline page.
