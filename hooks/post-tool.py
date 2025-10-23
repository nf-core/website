#!/usr/bin/env python3
"""
Claude Code Post-Tool Hook

Automatically runs after Claude uses Write, Edit, or MultiEdit tools.
Provides real-time formatting with prettier and prose quality checking with Vale.
"""

import sys
from pathlib import Path

# Add hooks directory to path for imports
hooks_dir = Path(__file__).parent
sys.path.insert(0, str(hooks_dir))

from utils import (
    read_hook_input, write_json_output, print_header, print_separator,
    get_file_paths_from_tool_input, format_file_with_prettier,
    lint_file_with_vale, should_process_file, should_lint_with_vale
)


def main():
    """Main hook execution."""
    # Read input from Claude Code
    input_data = read_hook_input()

    if not input_data:
        # No input available - exit silently
        sys.exit(0)

    # Extract relevant information
    tool_name = input_data.get('tool_name', 'Unknown')
    tool_input = input_data.get('tool_input', {})
    tool_response = input_data.get('tool_response', {})

    # Get file paths to process
    file_paths = get_file_paths_from_tool_input(tool_input, tool_name)

    if not file_paths:
        # No files to process - exit silently
        sys.exit(0)

    # Track processing results
    processed_files = []
    formatting_results = []
    vale_results = []
    has_prose_issues = False

    print_header("Post-Tool Formatting & Linting")
    print(f"🔧 Tool used: {tool_name}")
    print("📁 Processing files:")

    # Process each file
    for file_path in file_paths:
        print(f"   → {file_path}")

        # Skip if file doesn't exist
        if not Path(file_path).exists():
            print(f"   ⚠️  File not found: {file_path}")
            continue

        processed_files.append(file_path)

        # Format with prettier
        if should_process_file(file_path):
            print(f"🎨 Formatting {file_path} with prettier...")
            format_result = format_file_with_prettier(file_path)
            formatting_results.append({
                'file': file_path,
                'result': format_result
            })

            if format_result['success']:
                print(f"✅ {format_result['message']}")
            else:
                print(f"⚠️  {format_result['message']}")
        else:
            print(f"⏭️  Skipping formatting (unsupported file type)")

        # Lint with Vale
        if should_lint_with_vale(file_path):
            print(f"📝 Checking prose quality with Vale...")
            vale_result = lint_file_with_vale(file_path)
            vale_results.append({
                'file': file_path,
                'result': vale_result
            })

            if vale_result['success']:
                if vale_result['issues']:
                    has_prose_issues = True
                    print(f"📊 Vale found prose suggestions:")
                    # Indent the Vale output
                    for line in vale_result['output'].split('\n'):
                        if line.strip():
                            print(f"   {line}")
                    print("")
                else:
                    print(f"✅ Vale: No prose issues found")
            else:
                print(f"⚠️  Vale: {vale_result['message']}")
        else:
            print(f"⏭️  Skipping Vale check (not a markdown file)")

        print()  # Add spacing between files

    # Summary
    if processed_files:
        print("🎯 File processing complete!")

        if has_prose_issues:
            print()
            print("📚 Vale found some prose suggestions above. These help ensure:")
            print("   • Consistent terminology (nf-core, Nextflow, etc.)")
            print("   • Clear, readable documentation")
            print("   • Professional writing style")
            print()
            print("💡 Consider reviewing and addressing these suggestions for better documentation quality.")
        else:
            all_vale_passed = all(
                not result['result'].get('issues', False)
                for result in vale_results
                if result['result'].get('success', False)
            )
            if all_vale_passed and vale_results:
                print("✨ All prose checks passed - excellent writing!")

        # Provide additional context to Claude if there are issues
        if has_prose_issues:
            additional_context = []
            additional_context.append("Note: The files you just modified have some prose quality suggestions from Vale:")

            for result in vale_results:
                if result['result'].get('issues'):
                    additional_context.append(f"- {result['file']}: Found terminology and style suggestions")

            additional_context.append("These suggestions help maintain consistent nf-core documentation standards.")

            # Return structured output with additional context for Claude
            output = {
                "hookSpecificOutput": {
                    "hookEventName": "PostToolUse",
                    "additionalContext": "\n".join(additional_context)
                }
            }
            write_json_output(output)
        else:
            # Exit successfully without additional feedback
            pass

    else:
        print("ℹ️  No eligible files found to process")

    print_separator()
    print("✅ Hook execution complete")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n❌ Hook interrupted", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Hook error: {str(e)}", file=sys.stderr)
        sys.exit(1)