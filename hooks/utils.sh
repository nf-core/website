#!/usr/bin/env bash
# Utility functions for Claude Code hooks

# Check if a tool is available
check_tool() {
    local tool="$1"
    if ! command -v "$tool" &> /dev/null; then
        echo "⚠️  $tool not found - skipping"
        return 1
    fi
    return 0
}

# Check if file should be processed based on extension
should_process_file() {
    local file="$1"
    case "$file" in
        *.md|*.mdx|*.js|*.ts|*.astro|*.svelte|*.css|*.scss|*.json|*.yml|*.yaml)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Check if file should be linted with Vale
should_lint_with_vale() {
    local file="$1"
    case "$file" in
        *.md|*.mdx)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Format file with prettier if available
format_file() {
    local file="$1"

    if ! check_tool "prettier"; then
        return 1
    fi

    if ! should_process_file "$file"; then
        return 1
    fi

    echo "🎨 Formatting $file with prettier..."
    if prettier --write --log-level=error "$file" 2>/dev/null; then
        echo "✅ Formatted successfully"
        return 0
    else
        echo "❌ Prettier formatting failed"
        return 1
    fi
}

# Lint file with Vale if available
lint_file() {
    local file="$1"

    if ! check_tool "vale"; then
        return 1
    fi

    if ! should_lint_with_vale "$file"; then
        return 1
    fi

    echo "📝 Checking prose quality with Vale..."

    # Run Vale and capture output
    local vale_output
    vale_output=$(vale --config=.vale.ini "$file" 2>&1)
    local vale_exit_code=$?

    if [ $vale_exit_code -eq 0 ]; then
        echo "✅ Vale: No prose issues found"
        return 0
    else
        echo "📊 Vale found prose suggestions:"
        echo "$vale_output" | sed 's/^/   /'
        echo ""
        echo "💡 These suggestions help improve documentation clarity and consistency."
        return 1
    fi
}

# Get list of recently modified files (in last 5 minutes)
get_recent_files() {
    find . \( -name "*.md" -o -name "*.mdx" -o -name "*.js" -o -name "*.ts" -o -name "*.astro" -o -name "*.yml" -o -name "*.yaml" \) -newermt '5 minutes ago' 2>/dev/null | head -10
}

# Pretty print a separator
print_separator() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

# Print hook header
print_header() {
    local hook_name="$1"
    echo ""
    echo "🤖 Claude Code Hook: $hook_name"
    print_separator
}