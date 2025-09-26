#!/usr/bin/env bash
# Claude Code User Prompt Submit Hook
# Runs when user submits a prompt
# Provides quick prose quality check on recently modified documentation

# Source utilities
source "$(dirname "$0")/utils.sh"

# Get user prompt from arguments
USER_PROMPT="$1"

# Only run if the prompt is related to documentation/writing
# Check for common documentation keywords
if ! echo "$USER_PROMPT" | grep -qi -E "(document|write|edit|readme|guide|tutorial|prose|content|markdown|mdx|vale|style|writing|docs)"; then
    # Not a documentation-related prompt, exit silently
    exit 0
fi

# Get recently modified documentation files
RECENT_FILES=$(find . \( -name "*.md" -o -name "*.mdx" \) -newermt '10 minutes ago' 2>/dev/null | head -5)

# If no recent files, exit silently
if [ -z "$RECENT_FILES" ]; then
    exit 0
fi

print_header "Documentation Quality Check"

echo "🔍 Checking recently modified documentation files for prose quality..."
echo ""

issues_found=false
files_checked=0

while IFS= read -r file; do
    [ -z "$file" ] && continue
    [ ! -f "$file" ] && continue

    # Skip API reference files (they're ignored in Vale config anyway)
    if echo "$file" | grep -q "api_reference"; then
        continue
    fi

    # Skip old events (they're ignored in Vale config anyway)
    if echo "$file" | grep -qE "events/(2018|2019|2020|2021|2022|2023|2024)/"; then
        continue
    fi

    files_checked=$((files_checked + 1))
    echo "📄 Checking: $file"

    if ! check_tool "vale"; then
        echo "   ⚠️  Vale not available for prose checking"
        continue
    fi

    # Quick Vale check
    vale_output=$(vale --config=.vale.ini "$file" 2>&1)
    vale_exit_code=$?

    if [ $vale_exit_code -eq 0 ]; then
        echo "   ✅ Prose quality looks good"
    else
        issues_found=true
        echo "   📝 Found prose suggestions:"
        echo "$vale_output" | sed 's/^/      /'
    fi
    echo ""
done <<< "$RECENT_FILES"

# Summary for Claude
if [ $files_checked -eq 0 ]; then
    echo "ℹ️  No recent documentation files found to check"
elif [ "$issues_found" = true ]; then
    echo "📊 Summary: Found some prose suggestions in recent documentation."
    echo "💡 Consider these suggestions to improve clarity and consistency."
    echo "🎯 Common improvements: terminology consistency, readability, style."
else
    echo "✨ Summary: All recent documentation looks great!"
    echo "🎉 No prose issues found in recently modified files."
fi

print_separator
echo "🤖 Ready to help with your documentation request!"
echo ""