---
title: Bundled documentation
subtitle: Pipeline documentation must be hosted on the nf-core website
menu:
  main:
    weight: 60
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

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

Additional markdown pages (e.g., tutorials, FAQs) MAY be added under directories called `docs/usage/` or `docs/output/`, and will be automatically rendered on the nf-core website pipeline page as sub-pages of the corresponding section. These will be only listed in the sidebar. Providing links to the sub-pages within the main `usage.md` or `output.md` is highly recommended.
Markdown files elsewhere in `docs/` are not rendered on the website: such a page displays correctly on GitHub but returns a 404 on nf-co.re.
