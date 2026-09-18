---
title: CODEOWNERS
subtitle: Set up code ownership for your pipeline
---

A `CODEOWNERS` file assigns reviewers to pull requests based on which files changed.
GitHub reads the file, matches changed paths against patterns, and requests reviews from the listed owners.

Place the file at `.github/CODEOWNERS` in your pipeline repository.

:::note
CODEOWNERS is opt-in for nf-core pipelines.
It is not required by the nf-core template, but pipelines are free to use it.
:::

## Pattern matching

Each line is a file pattern followed by one or more GitHub usernames or team names.
Patterns follow most of the same rules as [gitignore](https://git-scm.com/docs/gitignore#_pattern_format) files.

- Order matters: the last matching pattern takes precedence.
- `*` matches any file but not `/`.
- `**` matches zero or more directories.
- `/` at the end matches only files directly in that directory.
- Paths are case sensitive.

:::warning
Some gitignore syntax rules do **not** work in CODEOWNERS:

- Escaping a pattern starting with `#` using `\`
- Using `!` to negate a pattern
- Using `[ ]` to define a character range

:::

## Common patterns in nf-core

Pipelines in nf-core use different levels of granularity.

### Minimal

A single catch-all rule that assigns all files to the pipeline maintainers.
This is the simplest and most common pattern:

```github-actions
# Maintainers
* @maintainer1 @maintainer2
```

Used by: [ampliseq](https://github.com/nf-core/ampliseq/blob/main/.github/CODEOWNERS), [mag](https://github.com/nf-core/mag/blob/main/.github/CODEOWNERS), [methylseq](https://github.com/nf-core/methylseq/blob/main/.github/CODEOWNERS), [pixelator](https://github.com/nf-core/pixelator/blob/main/.github/CODEOWNERS), [taxprofiler](https://github.com/nf-core/taxprofiler/blob/main/.github/CODEOWNERS).

### Granular

Ownership split by local modules and subworkflows, optionally down to individual tools with their own tests, configs, and test data.

```github-actions
# Maintainers
* @maintainer1 @maintainer2
*.nf.test* @maintainer1 @maintainer2
.github/workflows/ @nf-core/a-team

# Local modules
modules/local/* @maintainer1

# Local subworkflows
subworkflows/local/* @maintainer2
```

For pipelines with multiple tools maintained by different people, add per-tool rules:

```github-actions
# Tool A
modules/nf-core/tool-a/ @tool-a-maintainer
subworkflows/local/tool_a_workflow/ @tool-a-maintainer
tests/pipeline/default/tool_a*.nf.test @tool-a-maintainer
conf/test_tool_a.config @tool-a-maintainer

# Tool B
modules/nf-core/tool-b/ @tool-b-maintainer
subworkflows/local/tool_b_workflow/ @tool-b-maintainer
tests/pipeline/default/tool_b*.nf.test @tool-b-maintainer
conf/test_tool_b.config @tool-b-maintainer
```

Used by: [demultiplex](https://github.com/nf-core/demultiplex/blob/master/.github/CODEOWNERS), [rnafusion](https://github.com/nf-core/rnafusion/blob/master/.github/CODEOWNERS), [sarek](https://github.com/nf-core/sarek/blob/master/.github/CODEOWNERS).

## References

- [GitHub CODEOWNERS documentation](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-code-owners)
- [nf-core/proposals #40: CODEOWNER for pipelines](https://github.com/nf-core/proposals/issues/40)
