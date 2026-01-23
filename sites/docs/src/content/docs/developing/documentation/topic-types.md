---
title: Topic types
subtitle: Documentation topic types
shortTitle: Topic types
weight: 3
---

Documentation topic types provide a framework for organizing content.
Using consistent topic types helps readers find information quickly and ensures comprehensive coverage of features.

nf-core documentation follows the CTRT framework: Concept, Task, Reference, and Troubleshooting.
Most pages combine multiple topic types—for example, starting with a concept introduction followed by task instructions.

:::tip
Need help or have questions about documentation?
Contact the nf-core docs-team on Slack at **#team-docs**.
We're here to help!
:::

<!-- TODO: Add links to nf-core Slack -->

## Concept

Concept topics explain ideas, background information, and how things work.
They answer "what" something is and "why" you would use it.

Concept topics answer questions like:

- What is this feature?
- Why does it matter?
- When should I use it?
- How does it work at a high level?

Use concept topics to:

- Introduce new features or workflows
- Explain underlying principles or architecture
- Provide context before task instructions
- Describe relationships between components

### Writing concept topics

**Title format**: Use nouns rather than gerunds. For example, prefer "Object migration" over "Migrating objects" or "Migrate objects." Avoid generic titles like "Overview" or "Introduction"—use specific, searchable terms instead.

**Structure**: Keep concept topics focused on a single idea. Write one or two paragraphs that explain the concept and its value comprehensively for someone encountering it for the first time.

**What to avoid**: Don't include step-by-step instructions—save those for task topics. Avoid linking to procedures within concept content. If covering multiple related concepts, create separate topics rather than combining them.

## Task

Task topics provide step-by-step instructions for completing specific actions.
They guide users through procedures with clarity and consistency.

Task topics answer questions like:

- How do I accomplish this goal?
- What steps do I follow?
- What commands do I run?
- What options are available?

Use task topics to:

- Guide users through specific procedures
- Provide clear, actionable instructions
- Show command examples with expected output
- Explain prerequisites and assumptions

### Writing task topics

**Title format**: Begin with an active verb followed by a noun. For example, "Create a pipeline" or "Configure reference genomes." This structure clearly indicates the action required.

**Standard structure**:

1. **Introduction**: Brief explanation of when and why to perform this task
2. **Prerequisites** (optional): Required permissions, dependencies, or prior setup
3. **Steps**: Numbered list of actions to complete
4. **Results** (optional): Expected outcome and next steps

**Single-step tasks**: Format simple one-step procedures as bulleted list items rather than numbered steps. This visual distinction helps readers quickly identify straightforward tasks.

**Multiple approaches**: When documenting alternative methods, present the most common approach first. Use nested headings to organize different methods if all are necessary.

**What to avoid**: Don't include lengthy field-by-field explanations in task steps—link to reference content instead. Write in imperative mood (e.g., "Run the command" not "You can run the command").

## Reference

Reference topics provide detailed, structured information for quick lookup.
They present information in scannable formats like tables or lists, functioning like dictionary or encyclopedia entries.

Reference topics answer questions like:

- What parameters are available?
- What values are valid?
- What are the technical specifications?
- What configuration options exist?

Use reference topics to:

- Document all available parameters or options
- List valid values and defaults
- Provide complete API or configuration details
- Create tables of settings or specifications

### Writing reference topics

**Title format**: Use descriptive nouns that clearly identify the topic's subject matter. For example, "Pipeline parameters" or "Configuration options."

**Structure**: Begin with a brief opening sentence, then present content in organized tables or lists. Transform scattered details into cohesive, easy-to-navigate formats.

**Formatting priorities**: Prioritize scannability and searchability. Use tables, lists, and consistent formatting. Minimize explanatory text—reference topics prioritize completeness over narrative.

**What to avoid**: Avoid generic titles like "Important notes" or "Limitations." Instead, integrate limitation information contextually near related content, as these details often represent prerequisite information about functionality rather than standalone warnings.

## Troubleshooting

Troubleshooting topics address specific problems and their solutions.
They should appear as the final sections on a page.
When five or more troubleshooting topics exist, create a dedicated troubleshooting page instead.

Troubleshooting topics answer questions like:

- Why am I seeing this error?
- How do I fix this problem?
- What causes this issue?
- What are common mistakes?

Use troubleshooting topics to:

- Document known issues and solutions
- Explain error messages and their causes
- Provide diagnostic steps
- List common pitfalls and how to avoid them

### Writing troubleshooting topics

**Format types**: Troubleshooting content can take three forms:

1. **Introductory**: Introduces the troubleshooting section with context about potential issues
2. **Task-based**: Uses action-oriented titles like "Run debug tools" or "Verify syntax"
3. **Reference**: Presents specific errors with explanations and solutions

**Title format**: For error reference topics, include at least partial error messages in the title. Begin with severity indicators like "Error:" or "Warning:" when applicable. Include full error messages in the body text if the title is abbreviated.

**Solutions vs. workarounds**: Use "workaround" for temporary fixes and "resolution" or "resolve" for permanent solutions.

**Safety considerations**: When troubleshooting involves commands that modify data, include warnings about potential risks if commands aren't run correctly.

**Comprehensive coverage**: Include all troubleshooting information, even for rare or unlikely scenarios. The benefits of comprehensive coverage outweigh potential risks.

## Other topic types

Beyond the four primary types, additional specialized topic types serve specific documentation needs.

### Top-level pages

Top-level pages serve as entry points for major documentation sections.
They sit at the highest hierarchy level in navigation and introduce broad workflows.

**Purpose**: Provide navigation hubs that direct users to detailed content at the next level.

**Title format**: Use active verbs describing workflows, such as "Manage your infrastructure" or "Organize work with projects."

**Structure**:

- Brief workflow introduction
- Card-based listings of immediate child pages
- Links only to pages one level below

**What to include**: Focus on navigation and orientation rather than comprehensive documentation. Keep introductions brief and direct users to more detailed content.

### Get started pages

Get started pages introduce high-level concepts for broad feature areas.
They explain how multiple features work together within larger workflows.

**Purpose**: Provide conceptual understanding of how features interconnect, distinguishing from tutorials by focusing on concepts rather than task completion.

**When to use**: Only at the highest navigation levels for major workflow areas.

**Title format**: Use "Get started with [topic]" for page titles.

**Structure**:

1. Introduction explaining how features interconnect
2. Workflow diagram showing the process visually
3. Numbered steps organizing features by workflow
4. Brief descriptions with minimal inline links
5. "For more information" sections with relevant resources

**What to avoid**: Don't include detailed documentation—full details exist elsewhere. Avoid inline body links, reserve references for dedicated information sections.

### Tutorials

Tutorials provide end-to-end walkthroughs of complex workflows combining multiple tasks.
They guide users through complete processes spanning various features or tools.

**Purpose**: Help users achieve specific goals through sequential, multi-step procedures. Tutorials combine multiple tasks into comprehensive learning experiences.

**When to use**: When workflows involve complex, multi-step procedures that span various features. Not appropriate for single procedures—use task topics instead.

**Title format**: Begin with "Tutorial:" followed by an active verb, such as "Tutorial: Create a website" or "Tutorial: Deploy a pipeline."

**Structure**:

1. Introduction explaining outcomes
2. "Before you begin" section with prerequisites
3. Task sections with numbered instructions
4. Cross-references between related sections

**Writing style**: Use friendlier, more conversational language than standard documentation. Include encouraging phrases and use a supportive tone. Schematics are encouraged to clarify complex processes.

**Key requirement**: Provide a working example that readers can create themselves. Tutorials should result in tangible, functional outcomes.
