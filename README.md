<img src="public/images/logo/nf-core-logo.svg#gh-light-mode-only" width="400">
<img src="public/images/logo/nf-core-logo-darkbg.svg#gh-dark-mode-only" width="400">

# [nf-co.re](https://github.com/nf-core/website)

This repository contains code for the nf-core website: **<http://nf-co.re/>**

## Packages used

Here's how the website is built:

- Language: Javascript
- Frameworks:
  - [Astro](https://astro.build/) (static site generator),
  - [Svelte](https://svelte.dev/) (interactive components),
  - [Bootstrap 5.3](https://getbootstrap.com/docs/5.3/) (CSS framework)
- Tools:
  - [npm](https://www.npmjs.com/) (package manager)

## Development

### Getting the code

To make edits to the website, fork the repository to your own user on GitHub and then clone to your local system.

```bash
gh repo fork nf-core/website nf-core_website
cd nf-core_website/
```

### Installing dependencies

The website is built using [Astro](https://astro.build/), a static site generator.
To install the dependencies, run:

```bash
npm install
```

### Running a local server

Ok, you're ready! To run the website locally, just start astro dev mode:

```bash
npm run dev
```

You should then be able to access the website in your browser at [http://localhost:4321/](http://localhost:4321/).

### File structure

We follow Astro's [file structure](https://docs.astro.build/guides/project-structure) for the website.
The main files are:

- `src/pages/` - Astro pages
- `src/content/` - [Astro content collections](https://docs.astro.build/en/guides/content-collections/) (markdown files for events, docs, blog)
- `src/components/` - Astro/Svelte components
- `src/layouts/` - HTML layouts
- `src/styles/` - (S)CSS stylesheets
- `public/` - Static files (images, json files etc)

## Adding an event

To add an event, create a new markdown (or .mdx) file in `src/content/events/` with the following frontmatter:

```yaml
title: "Event Title"
subtitle: "A brief overview of the event"
type: "talk"  # Can be "talk", "hackathon", "training", "bytesize"
start_date: "YYYY-MM-DD"
end_date: "YYYY-MM-DD"
start_time: "HH:MM"
end_time: "HH:MM"
announcement:
  text: "Text on the announcement banner" # (optional)
  start: "YYYY-MM-DDTHH:MM:SS+HH:MM" # (required if announcement.text is used)
  end: "YYYY-MM-DDTHH:MM:SS+HH:MM" # (required if announcement.text is used)
location_name: "Name of the location" # (optional)
location_url: "URL to the location or to the section in the text with location description (e.g. `#gather-town`)" # (optional)
location_latlng: [48.2082, 16.3738] # Latitude and longitude of the location as an array " (optional)
address: "Address of the location" (optional)
duration: "Duration of the event in days" (optional)
embed_at: "in case this should be shown in the sidebar of a pipeline page (e.g. for a bytesize talk about the pipeline)" (optional)
import_typeform: true # If true, the event will be imported from a Typeform (see below)
```

## Adding a blog post

To add a blog post, create a new markdown (or mdx) file in `src/content/blog/` with the following frontmatter:

```yaml
title: "Your Blog Post Title"
subtitle: "A brief overview of your post's content"
headerImage: "Direct URL to an optional header image" (optional)
headerImageAlt: "Descriptive alt text for the header image (mandatory if a header image is used)"
pubDate: "Scheduled publication date and time (the post will go live post-website rebuild if the current date surpasses this timestamp). Format: YYYY-MM-DDTHH:MM:SS+HH:MM" (without quotes!)
authors: ["Author's Name"]  // Use a list format even if there is only one author.
label: ["Category1", "Category2"]  // This is optional and can include multiple categories.
announcement:
  text: "Text on the announcement banner" # (optional)
  start: "YYYY-MM-DDTHH:MM:SS+HH:MM" # (required if announcement.text is used)
  end: "YYYY-MM-DDTHH:MM:SS+HH:MM" # (required if announcement.text is used)
```

### Adding an announcement banner

You can show a short announcement banner on the website by adding additional information to the frontmatter of either a file inside `src/content/blog` or `src/content/events`. The following fields are available:

```yaml
announcement-text: "Your announcement text"
announcement-start: "Start date and time of the announcement. Format: YYYY-MM-DDTHH:MM:SS+HH:MM" (without quotes!)
announcement-end: "End date and time of the announcement. Format: YYYY-MM-DDTHH:MM:SS+HH:MM" (without quotes!) This is an optional field for events, where the start date of the event is the end date of the announcement by default.
```

### Updating the JSON files and cached markdown

Much of the site is powered by the JSON files in `/public` and the cached markdown files (from the pipeline docs) in `/.cache`.

They come pre-built with the repository, but if you want to rebuild them then you'll need to run the following commands. Note that you need to add a GITHUB_TOKEN inside a `.env` file to avoid hitting API limits (too early). See [instructions on how to get a GitHub OAuth token](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line) (the token only needs the `public_repo` permission).

```bash
npm run build-pipeline-json
npm run build-component-json
npm run build-cache-force
```

### Tools API docs

Tools docs are built using GitHub Actions on the nf-core/tools repo using Sphinx.
[These actions](https://github.com/nf-core/tools/blob/master/.github/workflows/tools-api-docs-release.yml) sync the built HTML files via FTP.

## Contribution guidelines

If you are looking forward to contribute to the website or add your institution to the official list of contributors, please have a look at the [CONTRIBUTING.md](./.github/CONTRIBUTING.md).

### Crafting a Blog Post

To publish a new blog post on the website, you'll need to create a Markdown file within the `src/content/blog/` directory. In this file, include the following frontmatter at the beginning:

```yaml
---
title: "Your Blog Post Title"
subtitle: "A brief overview of your post's content"
headerImage: "Direct URL to an optional header image"
headerImageAlt: "Descriptive alt text for the header image (mandatory if a header image is used)"
pubDate: "Scheduled publication date and time (the post will go live post-website rebuild if the current date surpasses this timestamp). Format: YYYY-MM-DDTHH:MM:SS.000+HH:MM"
authors: ["Author's Name"]  // Use a list format even if there is only one author.
label: ["Category1", "Category2"]  // This is optional and can include multiple categories.
---
```

> [!NOTE]
> The blog post will be visible on the website only if a rebuild of the site occurs after the date and time specified in the `pubDate` field.

By default the first paragraph of the blog post will be used as the preview text on the blog page. If you want to use a different paragraph, add the following comment after the paragraph you want to use:

```markdown
<!-- end of excerpt -->
```

or for MDX

<!-- prettier-ignore-start -->
```mdx
/* end of excerpt */
```
<!-- prettier-ignore-end -->

## Community

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## Credits

Phil Ewels ([@ewels](http://github.com/ewels/)) built the initial website, but there have been many contributors to the content and documentation.
Matthias Hörtenhuber ([@mashehu](https://github.com/mashehu)) worked on the concept and code for the new website rewrite.

See the [repo contributors](https://github.com/nf-core/website/graphs/contributors) for more details.

Kudos to the excellent [npm website](https://www.npmjs.com), which provided inspiration for the design of the pipeline pages.
