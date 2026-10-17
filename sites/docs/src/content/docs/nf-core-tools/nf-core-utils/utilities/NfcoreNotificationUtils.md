---
# Copied from docs/utilities/NfcoreNotificationUtils.md of nf-core/nf-core-utils 0.5.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "Pipeline Notifications Tutorial"
shortTitle: "Pipeline Notifications Tutorial"
weight: 6
---

# Pipeline Notifications Tutorial

Keeping users informed about pipeline execution is crucial for building user-friendly bioinformatics workflows. This tutorial teaches you how to implement comprehensive notification systems using the `NfcoreNotificationUtils` utility, covering everything from basic terminal summaries to sophisticated email and messaging integrations.

## 1. Overview

Modern bioinformatics pipelines often run for hours or days, making it essential to notify users about completion status, errors, and results. Users need to know:

- When their pipeline completes successfully or fails
- What parameters were used and what results were generated
- Where to find outputs and reports
- Performance metrics and resource usage

The `NfcoreNotificationUtils` utility provides three core notification channels:

- **Terminal Summaries**: Immediate, color-coded status information
- **Email Notifications**: Detailed reports with attachments and metadata
- **Instant Messaging**: Slack/Teams integration for team collaboration

### 1.1. Critical Session Context Requirement

:::warning{title="Session Context Required"}
**ALL notification functions require Nextflow session context and MUST be called from `workflow.onComplete` or `workflow.onError` handlers only.** They will fail with null pointer exceptions if called from the main workflow block.
:::

This requirement exists because notification functions need access to workflow metadata (success status, timing, resource usage) that is only available during completion handlers.

## 2. Getting Started with Notifications

Let's start with the simplest possible notification - a basic completion summary.

### 2.1. Your First Notification

Here's a minimal example that shows when and how to add notifications:

```nextflow title="basic_notification.nf"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import the completion summary function
include { completionSummary } from 'plugin/nf-core-utils'

params.input = "samples.csv"
params.outdir = "results"

workflow {
    log.info "Processing ${params.input}..."

    // Your pipeline logic here
    log.info "Analysis complete"
}

// CORRECT - Notification in completion handler
workflow.onComplete {
    completionSummary(params.monochrome_logs)
}
```

```console title="Terminal Output"
N E X T F L O W  ~  version 25.04.0
Launching `basic_notification.nf` [focused-darwin] - revision: abc1234

INFO [main] - Processing samples.csv...
INFO [main] - Analysis complete

-[nf-core/example] Pipeline completed successfully-
Completed at: 2024-01-15T10:30:45.123Z
Duration    : 45s
CPU hours   : 0.1 (95.2% cached)
Succeeded   : 3
```

### 2.2. Understanding the Session Context Problem

Let's understand why notification functions must be in completion handlers:

```nextflow title="session_context_demo.nf"
#!/usr/bin/env nextflow

include { completionSummary } from 'plugin/nf-core-utils'

workflow {
    log.info "Pipeline starting..."

    // WRONG - This will cause a NullPointerException
    // completionSummary(params.monochrome_logs)  // Session is null here!

    log.info "Pipeline logic completed"
}

workflow.onComplete {
    // CORRECT - Session context is available here
    log.info "Session available: ${session != null}"
    log.info "Workflow success: ${workflow.success}"
    log.info "Workflow duration: ${workflow.duration}"

    // Now the notification function works
    completionSummary(params.monochrome_logs)
}
```

### 2.3. Progressive Notification Building

Let's build up from basic terminal output to comprehensive notifications:

```nextflow title="progressive_notifications.nf" {15-23}
#!/usr/bin/env nextflow

include {
    completionSummary;
    completionEmail;
    imNotification
} from 'plugin/nf-core-utils'

params.email = null
params.hook_url = null

workflow {
    log.info "Starting comprehensive pipeline with notifications"
}

workflow.onComplete {
    // Level 1: Always show terminal summary
    completionSummary(params.monochrome_logs)

    // Level 2: Email notifications (if configured)
    if (params.email) {
        log.info "Sending completion email..."
        // We'll implement this next
    }

    // Level 3: Team notifications (if configured)
    if (params.hook_url) {
        log.info "Sending team notification..."
        // We'll implement this next
    }
}
```

## 3. Core Notification Functions

Now let's explore each notification function in detail, building from simple to complex.

### 3.1. completionSummary - Terminal Status Display

This function provides immediate visual feedback about pipeline completion status.

**Basic Terminal Summary**

The simplest notification shows essential completion information:

```nextflow title="terminal_summary_example.nf"
#!/usr/bin/env nextflow

include { completionSummary } from 'plugin/nf-core-utils'

workflow {
    log.info "Running analysis pipeline..."
    // Your pipeline processes here
}

workflow.onComplete {
    completionSummary(params.monochrome_logs)
}
```

**Function Parameters**

```groovy title="Function signature"
void completionSummary(boolean monochromeLogs = true)
```

| Parameter        | Type    | Default | Description                                           |
| ---------------- | ------- | ------- | ----------------------------------------------------- |
| `monochromeLogs` | Boolean | `true`  | If `true`, disables color codes for plain text output |

**Color-Coded Output Examples**

**Success with colors enabled:**

```console title="Successful completion (colors enabled)"
-[nf-core/rnaseq] Pipeline completed successfully-
Completed at: 2024-01-15T10:30:45.123Z
Duration    : 1h 23m 45s
CPU hours   : 45.2 (85.3% cached)
Succeeded   : 156
Cached      : 12
Ignored     : 0
```

**Failure with error details:**

```console title="Failed completion (colors enabled)"
-[nf-core/rnaseq] Pipeline completed with errors-
Completed at: 2024-01-15T09:15:23.456Z
Duration    : 15m 32s
CPU hours   : 2.1
Succeeded   : 23
Failed      : 1
Ignored     : 0

Exit status : 1
Error       : Process FASTQC (1) failed
```

**Advanced Terminal Summary Usage**

For pipelines that need custom completion reporting:

```nextflow title="advanced_terminal_summary.nf"
#!/usr/bin/env nextflow

include { completionSummary } from 'plugin/nf-core-utils'

workflow {
    // Pipeline logic...
}

workflow.onComplete {
    // Always show the summary
    completionSummary(params.monochrome_logs)

    // Add custom completion information
    if (workflow.success) {
        log.info """
        Analysis Summary:
        - Samples processed: ${workflow.stats.succeedCount}
        - Peak memory usage: ${workflow.stats.peakRss ?: 'N/A'}
        - Results saved to: ${params.outdir}
        """
    } else {
        log.error """
        Pipeline Failed:
        - Error: ${workflow.errorMessage ?: 'Unknown error'}
        - Failed tasks: ${workflow.stats.failedCount}
        - Check .nextflow.log for details
        """
    }
}
```

### 3.2. completionEmail - Detailed Email Reports

Email notifications provide comprehensive pipeline reports that users can refer to later.

**Understanding Email Structure**

Before implementing email notifications, let's understand what information needs to be organized:

```nextflow title="email_structure_demo.nf"
// The key to good emails is organizing parameters into logical groups
def summary_params = [
    'Core Nextflow options': [
        'runName': workflow.runName,
        'profile': workflow.profile,
        'revision': workflow.revision ?: 'N/A'
    ],
    'Input/output options': [
        'input': params.input ?: 'N/A',
        'outdir': params.outdir ?: 'N/A'
    ],
    'Resource usage': [
        'max_memory': params.max_memory,
        'max_cpus': params.max_cpus,
        'max_time': params.max_time
    ]
]
```

**Basic Email Implementation**

Let's implement a simple email notification:

```nextflow title="basic_email_notification.nf"
#!/usr/bin/env nextflow

include { completionEmail } from 'plugin/nf-core-utils'

params.email = null  // User provides email address
params.email_on_fail = null  // Optional separate failure email
params.outdir = "results"

workflow {
    log.info "Running pipeline that will send email notification..."
    // Your pipeline processes here
}

workflow.onComplete {
    // Only send email if user provided an email address
    if (params.email || params.email_on_fail) {

        // Organize parameters for the email
        def summary_params = [
            'Pipeline Information': [
                'runName': workflow.runName,
                'success': workflow.success,
                'duration': workflow.duration,
                'workDir': workflow.workDir
            ],
            'User Parameters': [
                'outdir': params.outdir,
                'profile': workflow.profile
            ]
        ]

        // Send the email
        completionEmail(
            summary_params,           // Parameter summary
            params.email,            // Success email address
            params.email_on_fail,    // Failure email address
            false,                   // Use HTML format (not plain text)
            params.outdir,           // Output directory
            params.monochrome_logs,  // Color settings
            []                       // MultiQC reports (empty for now)
        )

        log.info "Completion email sent"
    }
}
```

**Function Parameters Reference**

```groovy title="completionEmail signature"
void completionEmail(
    Map summaryParams,        // Grouped parameter summary
    String email,            // Primary email address
    String emailOnFail,      // Failure notification email
    boolean plaintextEmail,  // Use plain text instead of HTML
    String outdir,          // Output directory path
    boolean monochromeLogs, // Color settings for email content
    List multiqcReports    // MultiQC report attachments
)
```

**Advanced Email with MultiQC Reports**

For production pipelines, include MultiQC reports as attachments:

```nextflow title="advanced_email_notification.nf" focus={33-35,44}
#!/usr/bin/env nextflow

include { completionEmail } from 'plugin/nf-core-utils'

params.email = null
params.outdir = "results"

workflow {
    // Your pipeline processes that generate MultiQC reports
    MULTIQC(analysis_results)
}

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        def summary_params = [
            'Pipeline Summary': [
                'runName': workflow.runName,
                'success': workflow.success,
                'duration': workflow.duration,
                'exitStatus': workflow.exitStatus
            ],
            'Configuration': [
                'profile': workflow.profile,
                'container': workflow.containerEngine,
                'workDir': workflow.workDir
            ],
            'Results': [
                'outdir': params.outdir,
                'publishMode': params.publish_dir_mode
            ]
        ]

        // Include MultiQC reports as attachments
        def multiqc_report = file("${params.outdir}/multiqc/multiqc_report.html")
        def mqc_reports = multiqc_report.exists() ? [multiqc_report] : []

        completionEmail(
            summary_params,
            params.email,
            params.email_on_fail,
            params.plaintext_email ?: false,
            params.outdir,
            params.monochrome_logs,
            mqc_reports  // Attach MultiQC report if it exists
        )
    }
}
```

### 3.3. imNotification - Instant Messaging Integration

Team collaboration is enhanced with Slack and Microsoft Teams notifications that keep everyone informed about pipeline status.

**Understanding Webhook Integration**

Instant messaging notifications use webhooks to send JSON payloads to chat services:

```nextflow title="webhook_basics.nf"
// Slack webhook URL format
params.hook_url = "https://hooks.slack.com/services/YOUR/SLACK/WEBHOOK"

// Teams webhook URL format
// params.hook_url = "https://outlook.office.com/webhook/YOUR/TEAMS/WEBHOOK"
```

**Basic Slack/Teams Notification**

Let's implement basic instant messaging notifications:

```nextflow title="basic_im_notification.nf"
#!/usr/bin/env nextflow

include { imNotification } from 'plugin/nf-core-utils'

workflow.onComplete {
    // Show completion summary with colors (if enabled)
    completionSummary(params.monochrome_logs)
}

workflow.onError {
    // Show error summary
    completionSummary(params.monochrome_logs)
}
```

**Terminal Output Features:**

- **Color-coded Status**: Green for success, red for failure
- **Execution Statistics**: Duration, CPU time, memory usage
- **Process Summary**: Number of processes, cached tasks, ignored processes
- **Resource Information**: Peak memory, CPU efficiency
- **Configuration Summary**: Profile, container engine, work directory

**Example Output:**

```
-[nf-core/rnaseq] Pipeline completed successfully-
Completed at: 2024-01-15T10:30:45.123Z
Duration    : 1h 23m 45s
CPU hours   : 45.2 (85.3% cached)
Succeeded   : 156
Cached      : 12
Ignored     : 0
```

---

### `imNotification()`

**Description:**  
Sends JSON notifications to instant messaging webhooks (Slack, Microsoft Teams) with pipeline summary and workflow metadata. Supports rich formatting for better readability.

**Function Signature:**

```nextflow
void imNotification(Map summaryParams, String hookUrl)
```

**Parameters:**

- `summaryParams` (Map): Nested map of grouped summary parameters (same structure as `completionEmail`)
- `hookUrl` (String): Webhook URL for Slack, Teams, or other compatible service

**Usage Requirements:**

- Must be called from `workflow.onComplete` or `workflow.onError` handlers only
- Requires Nextflow session context (workflow metadata, timing, success status)
- Do not call from main workflow block
- Function handles null/empty hookUrl gracefully (logs warning and returns)

**Supported Webhook Types:**

- **Slack**: `https://hooks.slack.com/services/...`
- **Microsoft Teams**: `https://outlook.office.com/webhook/...`
- **Generic webhooks**: Any endpoint accepting JSON POST requests

**Usage Example:**

```nextflow
include { imNotification } from 'plugin/nf-core-utils'

workflow.onComplete {
    if (params.hook_url) {
        def summary_params = [
            'Pipeline Info': [
                'runName': workflow.runName,
                'success': workflow.success,
                'duration': workflow.duration,
                'exitStatus': workflow.exitStatus
            ],
            'Configuration': [
                'profile': workflow.profile,
                'container': workflow.containerEngine
            ]
        ]

        imNotification(summary_params, params.hook_url)
    }
}

workflow.onError {
    if (params.hook_url) {
        def summary_params = [
            'Error Details': [
                'runName': workflow.runName,
                'errorMessage': workflow.errorMessage,
                'exitStatus': workflow.exitStatus
            ]
        ]

        imNotification(summary_params, params.hook_url)
    }
}
```

**JSON Payload Structure:**

- **Service Detection**: Automatically detects Slack vs Teams format
- **Rich Formatting**: Uses service-specific formatting (markdown, cards)
- **Color Coding**: Green for success, red for failure
- **Comprehensive Data**: Includes all workflow metadata and custom parameters

---

## Complete Integration Example

Here's a comprehensive example showing all notification functions working together:

```nextflow expandable={62}
#!/usr/bin/env nextflow

include { completionSummary; completionEmail; imNotification } from 'plugin/nf-core-utils'

workflow {
    // Main workflow logic here
    log.info "Pipeline execution starting..."

    // Process data...

    log.info "Pipeline execution complete"
}

workflow.onComplete {
    // Always show terminal summary
    completionSummary(params.monochrome_logs)

    // Prepare summary parameters (reused for all notifications)
    def summary_params = [
        'Core Nextflow options': [
            'revision': workflow.revision ?: 'N/A',
            'runName': workflow.runName,
            'containerEngine': workflow.containerEngine,
            'profile': workflow.profile,
            'configFiles': workflow.configFiles.join(', ')
        ],
        'Input/output options': [
            'input': params.input ?: 'N/A',
            'outdir': params.outdir ?: 'N/A',
            'publish_dir_mode': params.publish_dir_mode
        ],
        'Reference genome options': [
            'genome': params.genome ?: 'None',
            'fasta': params.fasta ?: 'N/A'
        ],
        'Max job request options': [
            'max_cpus': params.max_cpus,
            'max_memory': params.max_memory,
            'max_time': params.max_time
        ]
    ]

    // Send completion email if configured
    if (params.email || params.email_on_fail) {
        completionEmail(
            summary_params,
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            multiqc_report.ifEmpty([])
        )
    }

    // Send instant message notification if configured
    if (params.hook_url) {
        imNotification(summary_params, params.hook_url)
    }

    log.info "All notifications sent successfully"
}

workflow.onError {
    // Show error summary
    completionSummary(params.monochrome_logs)

    // Send error notifications if configured
    if (params.email_on_fail || params.hook_url) {
        def error_params = [
            'Error Details': [
                'runName': workflow.runName,
                'errorMessage': workflow.errorMessage ?: 'Unknown error',
                'exitStatus': workflow.exitStatus,
                'errorReport': workflow.errorReport ?: 'No error report available'
            ],
            'Configuration': [
                'profile': workflow.profile,
                'workDir': workflow.workDir,
                'container': workflow.containerEngine
            ]
        ]

        // Send error email
        if (params.email_on_fail) {
            completionEmail(
                error_params,
                null,  // No regular email
                params.email_on_fail,  // Only send to error email
                params.plaintext_email,
                params.outdir,
                params.monochrome_logs,
                null  // No MultiQC report on error
            )
        }

        // Send error IM notification
        if (params.hook_url) {
            imNotification(error_params, params.hook_url)
        }
    }
}
```

---

## Best Practices

### 1. Parameter Organization

Group parameters logically for better email/notification formatting:

```nextflow
// GOOD - Grouped structure
def summary_params = [
    'Core Nextflow options': [
        'runName': workflow.runName,
        'profile': workflow.profile
    ],
    'Input/output options': [
        'input': params.input,
        'outdir': params.outdir
    ]
]

// AVOID - Flat structure
def summary_params = [
    runName: workflow.runName,
    profile: workflow.profile,
    input: params.input,
    outdir: params.outdir
]
```

### 2. Error Handling

Always provide fallback values for potentially null parameters:

```nextflow
def summary_params = [
    'Core Nextflow options': [
        'revision': workflow.revision ?: 'N/A',
        'runName': workflow.runName ?: 'Unknown',
        'profile': workflow.profile ?: 'default'
    ]
]
```

### 3. Testing Notifications

Test notification functions without sending actual notifications:

```nextflow
// Test mode - use mock webhook URL or disable actual sending
params.hook_url = params.test_mode ? 'https://httpbin.org/post' : params.hook_url
params.email = params.test_mode ? null : params.email
```

### 4. Security Considerations

- Store webhook URLs securely (use environment variables or secure config files)
- Don't log sensitive parameters in summary
- Use separate failure emails for security notifications

```nextflow
// Use environment variables for sensitive URLs
params.hook_url = System.getenv('SLACK_WEBHOOK_URL')

// Filter sensitive parameters from summary
def filtered_params = params.findAll { k, v ->
    !k.toLowerCase().contains('password') &&
    !k.toLowerCase().contains('token')
}
```

---

## Common Issues and Solutions

### Issue: "Session is null" Error

**Problem:**

```nextflow
// WRONG - This will cause null pointer exceptions
workflow {
    completionSummary(params.monochrome_logs)  // ERROR: session is null
}
```

**Solution:**

```nextflow
// CORRECT - Call from completion handlers
workflow.onComplete {
    completionSummary(params.monochrome_logs)  // Works correctly
}
```

### Issue: Email Not Sending

**Common Causes:**

1. SMTP configuration not set up in Nextflow config
2. Email parameters not provided
3. MultiQC report path incorrect

**Solution:**

```nextflow
// Check configuration
if (!params.email && !params.email_on_fail) {
    log.info "No email addresses provided - skipping email notification"
}

// Verify MultiQC report exists
def mqc_report = multiqc_report.ifEmpty([])
log.info "MultiQC reports to attach: ${mqc_report}"
```

### Issue: Webhook Failures

**Common Causes:**

1. Invalid webhook URL format
2. Network connectivity issues
3. Service-specific formatting problems

**Solution:**

```nextflow
// Validate webhook URL
if (params.hook_url && !params.hook_url.startsWith('http')) {
    log.warn "Invalid webhook URL format: ${params.hook_url}"
    params.hook_url = null
}

// Test webhook connectivity
if (params.hook_url) {
    log.info "Webhook configured: ${params.hook_url.take(20)}..."
}
```

---

## Testing

Test notification functions using nf-test:

```groovy
test("Notification functions work correctly") {
    when {
        params {
            email = "test@example.com"
            email_on_fail = "admin@example.com"
            hook_url = "https://hooks.slack.com/test"
            monochrome_logs = false
        }
    }

    then {
        assert workflow.success
        // Verify notification functions were called (check logs)
        assert workflow.stdout.contains("All notifications sent successfully")
    }
}

test("Error notifications work correctly") {
    when {
        // Force pipeline failure
        script """
        exit 1
        """
    }

    then {
        assert workflow.failed
        // Verify error notifications were sent
        assert workflow.stdout.contains("Error Details")
    }
}
```

---

## Integration with External Services

### Slack Integration

```nextflow
// Slack webhook configuration
params.hook_url = "https://hooks.slack.com/services/YOUR/SLACK/WEBHOOK"

// Slack-specific parameter formatting
def slack_summary = [
    'Pipeline Status': [
        'name': workflow.manifest.name,
        'version': workflow.manifest.version,
        'status': workflow.success ? 'SUCCESS' : 'FAILED'
    ]
]
```

### Microsoft Teams Integration

```nextflow
// Teams webhook configuration
params.hook_url = "https://outlook.office.com/webhook/YOUR/TEAMS/WEBHOOK"

// Teams supports rich card formatting
def teams_summary = [
    'Execution Summary': [
        'duration': workflow.duration,
        'cpu_hours': workflow.stats.computeTimeFmt,
        'memory_peak': workflow.stats.peakRss
    ]
]
```

For more information on setting up external service integrations, see the [nf-core documentation](https://nf-co.re/docs/usage/configuration#notifications).
