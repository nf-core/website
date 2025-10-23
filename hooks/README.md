# Claude Code Hooks for nf-core Website

This directory contains Claude Code hooks that automatically format and lint files when Claude modifies them, providing real-time feedback on prose quality using Python for robust JSON processing and integration.

## Available Hooks

### 🔧 `post-tool.py`

**Triggers**: After Claude uses Write, Edit, or MultiEdit tools
**Implementation**: Python 3 with structured JSON input/output
**Actions**:

- Reads Claude Code's structured JSON input via stdin
- Runs `prettier --write` on modified files for consistent formatting
- Runs `vale` prose linting on markdown/MDX files
- Provides immediate feedback on writing quality to Claude via JSON output
- Respects existing `.vale.ini` configuration
- Returns `additionalContext` to Claude when prose issues are found

### 📋 `user-prompt-submit.py`

**Triggers**: When user submits documentation-related prompts
**Implementation**: Python 3 with smart prompt detection
**Actions**:

- Reads user prompt from Claude Code's JSON input
- Checks recently modified documentation files (last 10 minutes)
- Runs Vale prose quality checks on recent files
- Only activates for documentation-related prompts (auto-detects keywords)
- Provides context to Claude via JSON output about documentation quality
- Adds guidance on nf-core documentation standards

### 🛠️ `utils.py`

**Purpose**: Shared Python utilities for all hooks
**Functions**:

- JSON input/output handling for Claude Code integration
- Tool availability checking (prettier, vale)
- File type filtering and validation
- Robust formatting and linting with error handling
- Recent file detection and processing
- Claude Code hook output formatting

## Setup Instructions

### 1. Configure Claude Code Settings

The hooks are automatically configured in `.claude/settings.local.json`:

```json
{
  "hooks": {
    "PostToolUse": [
      {
        "matcher": "Write|Edit|MultiEdit",
        "hooks": [
          {
            "type": "command",
            "command": "\"$CLAUDE_PROJECT_DIR\"/hooks/post-tool.py"
          }
        ]
      }
    ],
    "UserPromptSubmit": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "\"$CLAUDE_PROJECT_DIR\"/hooks/user-prompt-submit.py"
          }
        ]
      }
    ]
  }
}
```

### 2. Verify Dependencies

The hooks work best when these tools are installed:

- **Python 3** - Hook implementation language (system default)
- **prettier** - Code formatting (`npm install -g prettier`)
- **vale** - Prose linting (`brew install vale`)

### 3. Test the Setup

```bash
# Test post-tool hook with proper Claude Code JSON input
echo '{"tool_name":"Write","tool_input":{"file_path":"test.md"}}' | ./hooks/post-tool.py

# Test user-prompt-submit hook with documentation prompt
echo '{"prompt":"Help me write documentation"}' | ./hooks/user-prompt-submit.py

# Test user-prompt-submit hook with non-documentation prompt (should exit silently)
echo '{"prompt":"What is the weather?"}' | ./hooks/user-prompt-submit.py
```

## How It Works

### Post-Tool Hook Flow

1. **Claude modifies files** using Write/Edit/MultiEdit
2. **Hook triggers** automatically after tool execution
3. **Prettier formats** the modified files for consistency
4. **Vale checks prose** quality and terminology
5. **Feedback provided** to Claude about writing quality

### User-Prompt Hook Flow

1. **User submits prompt** containing documentation keywords
2. **Hook scans** for recently modified documentation files
3. **Vale checks** prose quality of recent changes
4. **Summary provided** before Claude responds

### Smart Filtering

- **File types**: Only processes `.md`, `.mdx`, `.js`, `.ts`, `.astro`, `.yml`, `.yaml`
- **Vale linting**: Only on markdown/MDX files
- **Respects `.vale.ini`**: Ignores API reference docs and old events
- **Performance**: Only processes changed files

## Benefits

### For Claude

- ✅ **Real-time feedback** on prose quality
- ✅ **Learning opportunity** from Vale suggestions
- ✅ **Consistent formatting** via prettier
- ✅ **Terminology awareness** (nf-core specific terms)

### For Contributors

- ✅ **Automatic formatting** of Claude-generated content
- ✅ **Prose quality assurance** built into the workflow
- ✅ **Consistent style** across all documentation
- ✅ **No additional setup** required

### For Project

- ✅ **Maintains high documentation standards**
- ✅ **Reduces review overhead** for maintainers
- ✅ **Ensures consistent terminology** usage
- ✅ **Professional documentation** quality

## Example Output

When Claude modifies a file with prose issues:

```
🤖 Claude Code Hook: Post-Tool Formatting & Linting
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🔧 Tool used: Edit
📁 Processing files:
   → docs/new-guide.md
🎨 Formatting docs/new-guide.md with prettier...
✅ Formatted successfully
📝 Checking prose quality with Vale...
📊 Vale found prose suggestions:
   7:8   error    Use 'nf-core' instead of 'nf_core'
   9:5   error    Use 'Bytesize' instead of 'bytesize'

💡 These suggestions help improve documentation clarity and consistency.
```

## Troubleshooting

### Hook Not Running

- Verify hooks are configured in Claude Code settings
- Check that hook files are executable (`chmod +x hooks/*.sh`)
- Ensure you're in the correct directory

### Missing Tools

- Install prettier: `npm install -g prettier`
- Install Vale: `brew install vale`
- Hooks will gracefully skip missing tools

### Vale Issues

- Hooks respect your existing `.vale.ini` configuration
- API reference docs and old events are automatically ignored
- Vale suggestions are informational, not blocking

This system creates a seamless feedback loop where Claude learns from its writing and continuously improves documentation quality! 🚀
