#!/usr/bin/env bash
# Claude Code Post-Tool Hook
# Runs after Claude uses Write, Edit, or MultiEdit tools
# Automatically formats files with prettier and checks prose with Vale

# Source utilities
source "$(dirname "$0")/utils.sh"

# Get hook input from Claude
HOOK_DATA="$1"

# Extract tool name and file paths from the hook data
TOOL_NAME=$(echo "$HOOK_DATA" | jq -r '.tool_name // empty' 2>/dev/null)
FILE_PATHS=$(echo "$HOOK_DATA" | jq -r '.file_paths[]? // empty' 2>/dev/null)

# If we can't parse JSON, try to extract file paths from environment or arguments
if [ -z "$FILE_PATHS" ] && [ -n "$CLAUDE_MODIFIED_FILES" ]; then
    FILE_PATHS="$CLAUDE_MODIFIED_FILES"
elif [ -z "$FILE_PATHS" ] && [ $# -gt 1 ]; then
    shift # Remove first argument (hook data)
    FILE_PATHS="$*"
fi

# If still no files, try to detect recently modified files
if [ -z "$FILE_PATHS" ]; then
    FILE_PATHS=$(get_recent_files)
fi

# Exit if no files to process
if [ -z "$FILE_PATHS" ]; then
    exit 0
fi

print_header "Post-Tool Formatting & Linting"

echo "🔧 Tool used: ${TOOL_NAME:-unknown}"
echo "📁 Processing files:"

# Convert FILE_PATHS to array and process each file
processed_any=false
vale_issues=false

# Handle both space-separated and newline-separated file lists
while IFS= read -r file; do
    # Skip empty lines
    [ -z "$file" ] && continue

    # Skip if file doesn't exist
    [ ! -f "$file" ] && continue

    echo "   → $file"

    # Format with prettier
    if format_file "$file"; then
        processed_any=true
    fi

    # Lint with Vale
    if ! lint_file "$file"; then
        vale_issues=true
    fi

    echo ""
done <<< "$(echo "$FILE_PATHS" | tr ' ' '\n')"

# Summary
if [ "$processed_any" = true ]; then
    echo "🎯 File processing complete!"

    if [ "$vale_issues" = true ]; then
        echo ""
        echo "📚 Vale found some prose suggestions above. These help ensure:"
        echo "   • Consistent terminology (nf-core, Nextflow, etc.)"
        echo "   • Clear, readable documentation"
        echo "   • Professional writing style"
        echo ""
        echo "💡 Consider reviewing and addressing these suggestions for better documentation quality."
    else
        echo "✨ All prose checks passed - excellent writing!"
    fi
else
    echo "ℹ️  No eligible files found to process"
fi

print_separator
echo "✅ Hook execution complete"