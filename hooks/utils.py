#!/usr/bin/env python3
"""
Utility functions for Claude Code hooks.

This module provides shared functionality for file processing, formatting,
and prose quality checking across all Claude Code hooks.
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
import shutil
import os
import time


def check_tool_available(tool: str) -> bool:
    """Check if a tool is available in the system PATH."""
    return shutil.which(tool) is not None


def should_process_file(file_path: str) -> bool:
    """Check if file should be processed based on extension."""
    file_path_obj = Path(file_path)
    return file_path_obj.suffix.lower() in {
        '.md', '.mdx', '.js', '.ts', '.tsx', '.jsx',
        '.astro', '.svelte', '.css', '.scss', '.json',
        '.yml', '.yaml'
    }


def should_lint_with_vale(file_path: str) -> bool:
    """Check if file should be linted with Vale."""
    file_path_obj = Path(file_path)
    return file_path_obj.suffix.lower() in {'.md', '.mdx'}


def format_file_with_prettier(file_path: str) -> Dict[str, Any]:
    """
    Format file with prettier if available.

    Returns:
        Dict with 'success' bool and 'message' str
    """
    if not check_tool_available('prettier'):
        return {
            'success': False,
            'message': 'prettier not available'
        }

    if not should_process_file(file_path):
        return {
            'success': False,
            'message': f'File type not supported: {Path(file_path).suffix}'
        }

    if not Path(file_path).exists():
        return {
            'success': False,
            'message': f'File does not exist: {file_path}'
        }

    try:
        result = subprocess.run(
            ['prettier', '--write', '--log-level=error', file_path],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            return {
                'success': True,
                'message': 'Formatted successfully'
            }
        else:
            return {
                'success': False,
                'message': f'Prettier failed: {result.stderr.strip()}'
            }
    except subprocess.TimeoutExpired:
        return {
            'success': False,
            'message': 'Prettier timed out'
        }
    except Exception as e:
        return {
            'success': False,
            'message': f'Prettier error: {str(e)}'
        }


def lint_file_with_vale(file_path: str) -> Dict[str, Any]:
    """
    Lint file with Vale if available.

    Returns:
        Dict with 'success' bool, 'issues' bool, 'output' str, and 'message' str
    """
    if not check_tool_available('vale'):
        return {
            'success': False,
            'issues': False,
            'output': '',
            'message': 'Vale not available'
        }

    if not should_lint_with_vale(file_path):
        return {
            'success': True,
            'issues': False,
            'output': '',
            'message': f'File type not supported for Vale: {Path(file_path).suffix}'
        }

    if not Path(file_path).exists():
        return {
            'success': False,
            'issues': False,
            'output': '',
            'message': f'File does not exist: {file_path}'
        }

    try:
        result = subprocess.run(
            ['vale', '--config=.vale.ini', file_path],
            capture_output=True,
            text=True,
            timeout=30
        )

        # Vale exit codes: 0 = no issues, >0 = issues found
        has_issues = result.returncode > 0

        return {
            'success': True,
            'issues': has_issues,
            'output': result.stdout.strip() if has_issues else '',
            'message': 'Vale check completed'
        }

    except subprocess.TimeoutExpired:
        return {
            'success': False,
            'issues': False,
            'output': '',
            'message': 'Vale timed out'
        }
    except Exception as e:
        return {
            'success': False,
            'issues': False,
            'output': '',
            'message': f'Vale error: {str(e)}'
        }


def get_recent_documentation_files(minutes: int = 10) -> List[str]:
    """
    Get recently modified documentation files.

    Args:
        minutes: How many minutes back to look for modifications

    Returns:
        List of file paths
    """
    cutoff_time = time.time() - (minutes * 60)
    recent_files = []

    try:
        # Find markdown and MDX files
        for pattern in ['**/*.md', '**/*.mdx']:
            for file_path in Path('.').glob(pattern):
                if not file_path.is_file():
                    continue

                # Skip API reference files (ignored in Vale config)
                if 'api_reference' in str(file_path):
                    continue

                # Skip old events (ignored in Vale config)
                if any(year in str(file_path) for year in [
                    'events/2018', 'events/2019', 'events/2020',
                    'events/2021', 'events/2022', 'events/2023', 'events/2024'
                ]):
                    continue

                # Check if file was modified recently
                if file_path.stat().st_mtime > cutoff_time:
                    recent_files.append(str(file_path))

        # Limit to 10 most recent files
        recent_files.sort(key=lambda f: Path(f).stat().st_mtime, reverse=True)
        return recent_files[:10]

    except Exception as e:
        print(f"Error finding recent files: {e}", file=sys.stderr)
        return []


def is_documentation_related_prompt(prompt: str) -> bool:
    """Check if prompt is related to documentation/writing."""
    doc_keywords = [
        'document', 'write', 'edit', 'readme', 'guide', 'tutorial',
        'prose', 'content', 'markdown', 'mdx', 'vale', 'style',
        'writing', 'docs', 'documentation'
    ]

    prompt_lower = prompt.lower()
    return any(keyword in prompt_lower for keyword in doc_keywords)


def print_header(title: str):
    """Print a formatted header."""
    print(f"\n🤖 Claude Code Hook: {title}")
    print("━" * 80)
    print()


def print_separator():
    """Print a separator line."""
    print("\n" + "━" * 80 + "\n")


def read_hook_input() -> Dict[str, Any]:
    """
    Read and parse JSON input from stdin.

    Returns:
        Parsed JSON data from Claude Code
    """
    try:
        if sys.stdin.isatty():
            # No input available
            return {}

        input_data = sys.stdin.read().strip()
        if not input_data:
            return {}

        return json.loads(input_data)

    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Error reading input: {e}", file=sys.stderr)
        return {}


def write_json_output(output: Dict[str, Any]):
    """Write JSON output to stdout for Claude Code."""
    try:
        print(json.dumps(output))
    except Exception as e:
        print(f"Error writing JSON output: {e}", file=sys.stderr)
        sys.exit(1)


def get_file_paths_from_tool_input(tool_input: Dict[str, Any], tool_name: str) -> List[str]:
    """
    Extract file paths from tool input based on tool type.

    Args:
        tool_input: The tool_input from Claude Code
        tool_name: The name of the tool that was used

    Returns:
        List of file paths to process
    """
    file_paths = []

    if tool_name in ['Write', 'Edit', 'MultiEdit']:
        # Single file operations
        if 'file_path' in tool_input:
            file_paths.append(tool_input['file_path'])
    elif tool_name == 'NotebookEdit':
        # Jupyter notebook operations
        if 'notebook_path' in tool_input:
            file_paths.append(tool_input['notebook_path'])

    # Filter to only existing files
    return [fp for fp in file_paths if Path(fp).exists()]