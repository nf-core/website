---
title: Style guide
subtitle: nf-core style guide rules
shortTitle: Style guide
weight: 3
---

<!-- TODO: Add high level section here -->

This guide covers the essential styling rules for nf-core documentation.

:::tip
Need help or have questions about documentation? Contact the nf-core docs-team on Slack at **#team-docs**. We're here to help!
:::

<!-- TODO: Add links to nf-core Slack -->

## Voice and tone

- Write in a **conversational but concise** style
- Use **active voice** and address readers directly ("you can configure..." not "the user can configure...")
- Focus on **user benefits**, not technical implementation details
- Avoid words like "simply", "easily", or "just" — they can be discouraging

## Grammar and punctuation

- Use **British English** spelling
- Write in **active voice** for clarity (exception: use passive voice when the subject is awkward)
- Include the **Oxford comma** in lists of three or more items
- Spell out **numbers zero through nine**; use numerals for 10 and above
- End complete sentences with a **period**
- Avoid **semicolons** — use separate sentences instead
- Avoid **em dashes** — use commas or separate sentences
- Keep sentences short and clear for international audiences

### Words to avoid

For clarity and translation:

- **"simply", "easily", "just"** — these can be discouraging
- **Ambiguous pronouns** like "it" — be specific if you can
- **Latin abbreviations** (e.g., i.e.) — spell them out ("for example", "that is")
- **"-ing" words** where possible — use direct verbs instead
- **"there is", "there are"** — these hide the subject
- **Idioms and cultural references** — these don't translate well

## Capitalisation

- Use **sentence case** for page titles and headings
- Match **UI text capitalisation** exactly as it appears
- Use **lowercase** for feature names unless they're proper nouns
- Follow the authoritative source for **third-party names** (GitHub, Docker, etc.)

## Text formatting

- **Bold** for UI elements: buttons, menu items, page names (match exact capitalisation)
- `Code formatting` for: filenames, paths, parameters, commands, variable names, error messages, user inputs
- Code blocks for multi-line examples and command-line instructions
- _Italics_ sparingly, primarily for emphasis

### Code blocks

- Use triple backticks with the language name for syntax highlighting
- Include blank lines before and after code blocks
- Use `<placeholder_text>` format for values users need to replace
- Example:

  ```bash
  nextflow run <pipeline_name> --input <input_file>
  ```

### Placeholder text

When showing examples with values users must replace:

- Use angle brackets: `<your-value>`
- Be descriptive: `<your_email@example.com>` not just `<email>`
- For test tokens in examples, use realistic fake values (see GitLab guide)

## Headings

- Use **sentence case** (not Title Case)
- Keep headings descriptive and scannable
- Don't skip heading levels (don't jump from H2 to H4)
- Avoid using H1 (`#`) in Markdown files — it's generated from the title
- Keep heading levels between **H2 and H4** — if you need H5, consider splitting into a new page

## Lists

- Keep list items **parallel in structure** (all start with verbs, or all are noun phrases)
- Use **numbered lists** for sequential steps
- Use **bulleted lists** for non-sequential items
- End list items with periods only if they're complete sentences

## Links

- Use **descriptive link text** — avoid "click here" or bare URLs
  - Good: `See the [pipeline parameters](link)`
  - Bad: `Click [here](link) for parameters`
- Don't duplicate the same link multiple times on one page
- Avoid links in headings
- Keep the entire link on one line (don't break across lines)
- For internal links, use relative paths
- For external links, use full URLs
- When linking to code, link to specific commits (not branches) so line numbers don't shift

## Tables

- Use tables to organise complex information
- **No empty cells** — use "N/A" or "None" instead
- Use **sentence case** for headers
- Keep columns consistently spaced
- Place description columns on the right when possible

## Images and screenshots

- Use sparingly — text is easier to maintain and translate
- Compress images to reduce file size
- Keep dimensions reasonable (max 1000px wide, 500px tall)
- Remove any personally identifiable information
- Use realistic examples (not placeholder text like "lorem ipsum")

### Alt text

- Keep it under 155 characters
- Use sentence case and end with a period
- Describe the **context**, not just the content
- Don't start with "Image of" or "Screenshot of"
- No formatting (bold, italics, or code) in alt text

## Notes and alerts

Use alert boxes sparingly for important information:

- **Note**: Additional context
- **Tip**: Helpful tips
- **Warning**: Potential issues or breaking changes
- **Important**: Critical information users must know

Don't stack multiple alerts consecutively — combine or restructure the content.

## Acronyms and abbreviations

- Spell out acronyms on first use per page: "Continuous Integration (CI)"
- Don't repeat the spelled-out version on the same page
- Avoid acronyms in page titles unless widely recognised
- Pluralise without apostrophes: "APIs" not "API's"

## Writing style tips

- **Get to the point** — don't write "This page explains..." — just explain it
- **Focus on facts** — avoid marketing language and subjective claims
- **Be specific** — instead of "improves performance", say "reduces build time by 50%"
- **Stay current** — don't promise future features; link to issues if relevant
- **Write for translation** — even if not translating, clear writing helps everyone
- **Break up long pages** — if a page gets too long, consider splitting into multiple pages
