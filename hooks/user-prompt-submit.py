#!/usr/bin/env python3
"""
Claude Code User Prompt Submit Hook

Runs when user submits a prompt. Provides quick prose quality check
on recently modified documentation and adds relevant context for Claude.
"""

import sys
from pathlib import Path

# Add hooks directory to path for imports
hooks_dir = Path(__file__).parent
sys.path.insert(0, str(hooks_dir))

from utils import (
    read_hook_input, write_json_output, print_header, print_separator,
    get_recent_documentation_files, lint_file_with_vale,
    is_documentation_related_prompt
)


def main():
    """Main hook execution."""
    # Read input from Claude Code
    input_data = read_hook_input()

    if not input_data:
        # No input available - exit silently
        sys.exit(0)

    # Extract the user prompt
    user_prompt = input_data.get('prompt', '')

    if not user_prompt:
        # No prompt provided - exit silently
        sys.exit(0)

    # Only run if the prompt is related to documentation/writing
    if not is_documentation_related_prompt(user_prompt):
        # Not a documentation-related prompt, exit silently
        sys.exit(0)

    # Get recently modified documentation files
    recent_files = get_recent_documentation_files(minutes=10)

    if not recent_files:
        # No recent documentation files, but still provide context about documentation intent
        additional_context = (
            "The user's prompt appears to be documentation-related. "
            "Consider documentation best practices: consistent terminology (nf-core, Nextflow), "
            "clear structure, and professional writing style."
        )

        output = {
            "hookSpecificOutput": {
                "hookEventName": "UserPromptSubmit",
                "additionalContext": additional_context
            }
        }
        write_json_output(output)
        sys.exit(0)

    print_header("Documentation Quality Check")
    print("🔍 Checking recently modified documentation files for prose quality...")
    print()

    issues_found = False
    files_checked = 0
    issue_summary = []

    for file_path in recent_files:
        if not Path(file_path).exists():
            continue

        files_checked += 1
        print(f"📄 Checking: {file_path}")

        # Quick Vale check
        vale_result = lint_file_with_vale(file_path)

        if vale_result['success']:
            if vale_result['issues']:
                issues_found = True
                print("   📝 Found prose suggestions:")
                # Show the Vale output with proper indentation
                for line in vale_result['output'].split('\n'):
                    if line.strip():
                        print(f"      {line}")
                issue_summary.append(f"• {file_path}: Found prose suggestions")
            else:
                print("   ✅ Prose quality looks good")
        else:
            print(f"   ⚠️  Vale not available: {vale_result['message']}")

        print()

    # Summary for user display
    if files_checked == 0:
        print("ℹ️  No recent documentation files found to check")
    elif issues_found:
        print("📊 Summary: Found some prose suggestions in recent documentation.")
        print("💡 Consider these suggestions to improve clarity and consistency.")
        print("🎯 Common improvements: terminology consistency, readability, style.")
    else:
        print("✨ Summary: All recent documentation looks great!")
        print("🎉 No prose issues found in recently modified files.")

    print_separator()
    print("🤖 Ready to help with your documentation request!")
    print()

    # Prepare additional context for Claude
    context_parts = []

    context_parts.append("Documentation context:")

    if issues_found:
        context_parts.append(f"Recent prose quality check found suggestions in {len([s for s in issue_summary])} files:")
        context_parts.extend(issue_summary)
        context_parts.append("Please consider nf-core documentation standards: use 'nf-core' (not 'nf_core'), 'Nextflow' (not 'nextflow'), proper capitalization, and clear writing.")
    else:
        context_parts.append(f"Recent prose quality check: all {files_checked} recently modified documentation files look good.")

    context_parts.append("Focus on maintaining consistent terminology and professional writing style.")

    # Return structured output with additional context
    output = {
        "hookSpecificOutput": {
            "hookEventName": "UserPromptSubmit",
            "additionalContext": "\n".join(context_parts)
        }
    }
    write_json_output(output)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n❌ Hook interrupted", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Hook error: {str(e)}", file=sys.stderr)
        sys.exit(1)