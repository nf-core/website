#!/usr/bin/env python3
"""
Split multi-sentence paragraphs in markdown files.

This script finds paragraphs with multiple sentences (ending in ., ?, or !)
and splits them onto consecutive lines so markdown still treats them as paragraphs.
"""

import os
import re
from pathlib import Path


def is_in_code_block(lines, line_idx):
    """Check if a line is inside a code block."""
    code_block_count = 0
    for i in range(line_idx):
        if lines[i].strip().startswith("```") or lines[i].strip().startswith("~~~"):
            code_block_count += 1
    return code_block_count % 2 == 1


def is_frontmatter(lines, line_idx):
    """Check if a line is in frontmatter (between --- markers at start of file)."""
    if line_idx == 0 and lines[0].strip() == "---":
        return False  # First line of frontmatter

    # Check if we're before the closing ---
    found_opening = False
    for i in range(line_idx):
        if lines[i].strip() == "---":
            if not found_opening:
                found_opening = True
            else:
                return False  # We're past the closing frontmatter

    return found_opening


def split_sentences_in_line(line):
    """
    Split a line with multiple sentences onto separate lines.
    Sentences end with ., ?, or ! (possibly followed by quotes, closing parens, etc.)
    """
    # Don't process empty lines or lines that are clearly markdown syntax
    if not line.strip() or line.strip().startswith(('#', '-', '*', '>', '|', '```', ':::')):
        return [line]

    # Don't process lines that look like numbered lists (e.g., "1. Something")
    if re.match(r'^\s*\d+\.\s', line):
        return [line]

    # Pattern to match sentence endings: . ? ! followed by optional quotes/parens and then a space
    # Use negative lookbehind to avoid matching numbered lists
    sentence_pattern = r'(?<!\d)([.?!]["\')\]]*)\s+(?=[A-Z"\'\(])'

    # Find all sentence boundaries
    parts = re.split(sentence_pattern, line)

    # Reconstruct sentences
    sentences = []
    current = ""
    for i, part in enumerate(parts):
        current += part
        # If this part is a sentence ending (. ? !) and next part starts sentence
        if i < len(parts) - 1 and re.match(r'^[.?!]["\')\]]*$', part.strip()):
            sentences.append(current.rstrip())
            current = ""

    # Add any remaining text
    if current.strip():
        sentences.append(current.rstrip())

    # Only return multiple lines if we found multiple sentences
    if len(sentences) > 1:
        return sentences
    return [line]


def process_file(filepath):
    """Process a single markdown file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    modified = False
    new_lines = []

    i = 0
    while i < len(lines):
        line = lines[i]

        # Skip if in code block or frontmatter
        if is_in_code_block(lines[:i+1], i) or is_frontmatter(lines, i):
            new_lines.append(line)
            i += 1
            continue

        # Try to split sentences
        split_lines = split_sentences_in_line(line.rstrip('\n'))

        if len(split_lines) > 1:
            # We found multiple sentences - add them as separate lines
            for split_line in split_lines:
                new_lines.append(split_line + '\n')
            modified = True
        else:
            new_lines.append(line)

        i += 1

    # Write back if modified
    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        return True
    return False


def main():
    """Find and process all markdown files."""
    docs_dir = Path.cwd()

    # Find all .md and .mdx files, excluding getting_started
    md_files = []
    for pattern in ['**/*.md', '**/*.mdx']:
        for filepath in docs_dir.glob(pattern):
            if 'getting_started' not in str(filepath) and 'split_sentences.py' not in str(filepath):
                md_files.append(filepath)

    print(f"Found {len(md_files)} markdown files to process")

    modified_count = 0
    for filepath in sorted(md_files):
        try:
            if process_file(filepath):
                modified_count += 1
                print(f"✓ Modified: {filepath.relative_to(docs_dir)}")
        except Exception as e:
            print(f"✗ Error processing {filepath.relative_to(docs_dir)}: {e}")

    print(f"\nProcessed {len(md_files)} files, modified {modified_count}")


if __name__ == "__main__":
    main()
