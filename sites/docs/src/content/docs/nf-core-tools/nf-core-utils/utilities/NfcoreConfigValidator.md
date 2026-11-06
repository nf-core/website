---
# Copied from docs/utilities/NfcoreConfigValidator.md of nf-core/nf-core-utils 0.7.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "Configuration Validation Tutorial"
shortTitle: "Configuration Validation Tutorial"
weight: 5
---

# Configuration Validation Tutorial

Configuration validation is the first line of defense against common pipeline execution issues. This guide teaches you how to implement robust configuration checking in your nf-core pipelines using the `NfcoreConfigValidator` utility.

## 1. Overview

Poor configuration is one of the leading causes of pipeline failures and user frustration. Users often:

- Run pipelines without custom profiles, leading to suboptimal resource allocation
- Make syntax errors in profile specifications
- Provide positional arguments that are silently ignored
- Use default settings that aren't suitable for their environment

The `NfcoreConfigValidator` utility addresses these issues by providing proactive validation that catches problems early and guides users toward best practices.

### 1.1. What You'll Learn

By the end of this tutorial, you'll understand how to:

- Implement configuration validation that prevents common user errors
- Provide helpful, color-coded feedback messages
- Guide users toward nf-core configuration best practices
- Handle edge cases gracefully without breaking pipelines

### 1.2. Core Validation Functions

The utility provides two essential validation functions:

- **`checkConfigProvided()`**: Validates that custom configuration or profiles are being used
- **`checkProfileProvided(args, monochromeLogs)`**: Validates command-line profile arguments and warns about common mistakes

These checks are explicit: the plugin's observer does not run them. Call each once
during pipeline initialization, or retain the equivalent template-local check.
Do not call both implementations of the same check. Pipelines that only load the
plugin must add the calls if they want validation, including trailing-comma
rejection. An uncaught profile validation error aborts Nextflow with a nonzero exit
status and a single `ERROR` label (with optional colour); missing custom config and
positional arguments remain warnings. Call validation from pipeline initialization,
not an observer callback, where Nextflow catches and logs exceptions.

## 2. Getting Started with Configuration Validation

Let's start with a simple example that demonstrates why configuration validation matters.

### 2.1. The Problem: Default Configuration Issues

Consider what happens when users run pipelines without proper configuration:

```bash title="Problematic user commands"
# User runs with default settings - resource waste and poor performance
nextflow run nf-core/rnaseq --input samples.csv --outdir results

# User provides conflicting arguments
nextflow run nf-core/rnaseq --input samples.csv -profile docker,

# User provides positional arguments that are ignored
nextflow run nf-core/rnaseq samples.csv --outdir results
```

Without validation, these commands might run but produce suboptimal results or unexpected behavior.

### 2.2. The Solution: Proactive Validation

Let's implement basic validation that catches these issues:

```nextflow title="basic_validation.nf"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import validation functions
include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'

// Validate configuration early - before any heavy processing
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

workflow {
    log.info "Configuration validation passed - proceeding with analysis"
    // Your pipeline logic here
}
```

```console title="Output with good configuration"
$ nextflow run basic_validation.nf -profile docker --input samples.csv

N E X T F L O W  ~  version 25.04.0
Launching `basic_validation.nf` [peaceful-euler] - revision: abc1234

INFO [main] - Configuration validation passed - proceeding with analysis
```

```console title="Output with configuration warnings"
$ nextflow run basic_validation.nf --input samples.csv

N E X T F L O W  ~  version 25.04.0
Launching `basic_validation.nf` [focused-darwin] - revision: abc1234

WARN [main] - No custom Nextflow configuration detected! Consider using an institutional config or custom profile.
INFO [main] - Configuration validation passed - proceeding with analysis
```

## 3. Core Validation Functions

### 3.1. checkConfigProvided - Custom Configuration Detection

This function encourages users to provide custom configurations for optimal performance.

**Understanding Configuration Types**

Nextflow can be configured in several ways, listed in order of preference:

1. **Institutional configs** (best): Pre-configured settings for specific institutions
2. **Custom profiles**: User-defined profiles for specific environments
3. **Custom config files**: Project-specific configuration files
4. **Default settings** (suboptimal): Basic Nextflow defaults

**Basic Implementation**

```nextflow title="config_detection_example.nf"
#!/usr/bin/env nextflow

include { checkConfigProvided } from 'plugin/nf-core-utils'

// Simple validation
def hasCustomConfig = checkConfigProvided()

if (hasCustomConfig) {
    log.info "✓ Custom configuration detected"
} else {
    log.warn "⚠ Using default configuration - consider adding a profile"
}

workflow {
    log.info "Starting pipeline with ${hasCustomConfig ? 'custom' : 'default'} configuration"
}
```

**Function Reference**

```groovy title="Function signature"
boolean checkConfigProvided()
```

| Aspect                   | Details                                             |
| ------------------------ | --------------------------------------------------- |
| **Parameters**           | None (uses session context automatically)           |
| **Returns**              | `boolean` - `true` if custom configuration detected |
| **Session Dependencies** | Requires Nextflow session context                   |
| **Side Effects**         | Logs warning messages if no custom config found     |

**Advanced Configuration Guidance**

For production pipelines, provide more detailed guidance:

```nextflow title="advanced_config_guidance.nf"
#!/usr/bin/env nextflow

include { checkConfigProvided } from 'plugin/nf-core-utils'

def hasCustomConfig = checkConfigProvided()

if (!hasCustomConfig) {
    log.warn """
    ==========================================
    NO CUSTOM CONFIGURATION DETECTED
    ==========================================

    For optimal performance, consider:

    1. Institutional configs:
       nextflow run pipeline.nf -profile myInstitution

    2. Container profiles:
       nextflow run pipeline.nf -profile docker
       nextflow run pipeline.nf -profile singularity

    3. Custom config files:
       nextflow run pipeline.nf -c my_config.config

    4. HPC profiles:
       nextflow run pipeline.nf -profile test,slurm,docker

    See: https://nf-co.re/docs/usage/configuration
    ==========================================
    """
}

workflow {
    // Pipeline continues regardless of configuration status
    log.info "Pipeline starting..."
}
```

### 3.2. checkProfileProvided - Profile Argument Validation

This function validates command-line profile arguments and catches common user mistakes before they cause issues.

**Understanding Common Profile Mistakes**

Users frequently make syntax errors when specifying profiles:

```bash title="Common profile mistakes"
# Trailing comma (causes parsing errors)
nextflow run pipeline.nf -profile docker,

# Space-separated values (only first is used)
nextflow run pipeline.nf -profile test docker

# Positional arguments (silently ignored)
nextflow run pipeline.nf samples.csv --outdir results

# Correct usage
nextflow run pipeline.nf -profile test,docker --input samples.csv --outdir results
```

**Basic Profile Validation**

Let's implement profile validation with helpful error messages:

```nextflow title="profile_validation_example.nf"
#!/usr/bin/env nextflow

include { checkProfileProvided } from 'plugin/nf-core-utils'

// Validate profile arguments - catches syntax errors early
checkProfileProvided(args, params.monochrome_logs)

workflow {
    log.info "Profile validation passed"
    log.info "Using profile: ${workflow.profile}"
}
```

**Color-Coded Error Messages**

The function supports color formatting for better user experience:

```nextflow title="colored_validation.nf" hl_lines="5"
#!/usr/bin/env nextflow

include { checkProfileProvided } from 'plugin/nf-core-utils'

// Enable colors for better error visibility (when supported)
checkProfileProvided(args, params.monochrome_logs)

workflow {
    log.info "Starting pipeline with profile: ${workflow.profile}"
}
```

**Without colors (default):**

```console title="Monochrome error output"
ERROR ~ The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!
HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.
```

**With colors enabled:**

```console title="Color error output"
[31mERROR[0m ~ The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!
[33mHINT[0m: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.
```

**Function Reference**

```groovy title="Function signature"
void checkProfileProvided(List args, boolean monochromeLogs = true)
```

| Parameter        | Type    | Default  | Description                                           |
| ---------------- | ------- | -------- | ----------------------------------------------------- |
| `args`           | List    | Required | Command-line arguments (use built-in `args` variable) |
| `monochromeLogs` | Boolean | `true`   | If `true`, disables color codes in error messages     |

**Comprehensive Profile Validation**

For production pipelines, implement comprehensive validation:

```nextflow title="comprehensive_profile_validation.nf"
#!/usr/bin/env nextflow

include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'

// Validate both configuration and profile arguments
def hasCustomConfig = checkConfigProvided()

if (!hasCustomConfig) {
    log.info """
    Consider using a configuration profile for better performance:
    - For Docker: -profile docker
    - For Singularity: -profile singularity
    - For HPC: -profile test,slurm,docker
    - For testing: -profile test
    """
}

// Validate profile syntax and warn about common mistakes
checkProfileProvided(args, params.monochrome_logs)

workflow {
    log.info """
    Pipeline Configuration Summary:
    - Custom config: ${hasCustomConfig ? 'Yes' : 'No (using defaults)'}
    - Profile: ${workflow.profile ?: 'none'}
    - Container engine: ${workflow.containerEngine ?: 'none'}
    """
}
```

## 4. Real-World Integration Patterns

Let's explore how to integrate configuration validation into different types of pipelines.

### 4.1. Basic nf-core Pipeline Integration

Here's the standard integration pattern for nf-core pipelines:

```nextflow title="standard_nfcore_pipeline.nf" hl_lines="5-8"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Essential validation at the top of main.nf
include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

// Pipeline parameters
params.input = null
params.outdir = "./results"
params.genome = null

// Import additional pipeline components
include { INPUT_CHECK } from './subworkflows/local/input_check'
include { FASTQC } from './modules/nf-core/fastqc/main'

workflow {
    if (!params.input) {
        error "Please provide an input sample sheet with --input"
    }

    log.info "Starting pipeline with validated configuration"

    // Pipeline workflow logic
    ch_input = INPUT_CHECK(file(params.input, checkIfExists: true))
    FASTQC(ch_input)
}
```

### 4.2. Advanced Validation with Environment Detection

For pipelines that need to adapt to different environments:

```nextflow title="environment_aware_validation.nf"
#!/usr/bin/env nextflow

include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'

// Detect environment and adjust validation accordingly
def isCI = System.getenv('CI') != null
def isHPC = System.getenv('SLURM_JOB_ID') != null
def isCloud = workflow.workDir.toString().startsWith('s3://')

// Validate configuration with environment-specific guidance
def hasCustomConfig = checkConfigProvided()
if (!hasCustomConfig) {
    if (isHPC) {
        log.warn "Running on HPC without custom config. Consider: -profile test,slurm,singularity"
    } else if (isCloud) {
        log.warn "Running on cloud without custom config. Consider: -profile test,aws,docker"
    } else {
        log.warn "Consider using -profile docker or -profile singularity for local execution"
    }
}

// Enable colors based on environment
checkProfileProvided(args, isCI ? true : params.monochrome_logs)

workflow {
    log.info """
    Environment Detection:
    - CI/CD: ${isCI}
    - HPC: ${isHPC}
    - Cloud: ${isCloud}
    - Config: ${hasCustomConfig ? 'Custom' : 'Default'}
    """
}
```

### 4.3. Validation with Error Recovery

Sometimes you want to provide suggestions and continue rather than fail:

```nextflow title="graceful_validation.nf"
#!/usr/bin/env nextflow

include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'

// Graceful validation that doesn't stop the pipeline
def performGracefulValidation() {
    try {
        def hasCustomConfig = checkConfigProvided()
        checkProfileProvided(args, params.monochrome_logs)

        if (!hasCustomConfig && !workflow.profile) {
            log.warn """
            ⚠️  PERFORMANCE WARNING ⚠️

            No custom configuration detected. This may lead to:
            - Inefficient resource usage
            - Slower execution times
            - Compatibility issues

            Recommended action: Add -profile docker or -profile singularity

            Pipeline will continue with default settings...
            """

            // Add a delay to ensure user sees the warning
            Thread.sleep(3000)
        }

        return true

    } catch (Exception e) {
        log.error "Validation error: ${e.message}"
        log.warn "Pipeline will attempt to continue, but may fail..."
        return false
    }
}

def validationPassed = performGracefulValidation()

workflow {
    log.info "Validation status: ${validationPassed ? 'Passed' : 'Warning'}"

    // Pipeline continues regardless of validation results
    log.info "Starting pipeline execution..."
}
```

## 5. Troubleshooting Common Issues

### 5.1. Profile Validation Errors

| Error                                 | Cause                                | Solution                                     |
| ------------------------------------- | ------------------------------------ | -------------------------------------------- |
| "cannot end with a trailing comma"    | `nextflow run -profile docker,`      | Remove trailing comma: `-profile docker`     |
| "positional argument detected"        | `nextflow run pipeline.nf input.csv` | Use proper syntax: `--input input.csv`       |
| "multiple values separated by spaces" | `-profile test docker`               | Use comma separation: `-profile test,docker` |

### 5.2. Configuration Detection Issues

**Problem**: `checkConfigProvided()` returns `false` even with custom config

**Debugging steps**:

```nextflow title="config_debug.nf"
include { checkConfigProvided } from 'plugin/nf-core-utils'

def hasConfig = checkConfigProvided()

// Debug configuration detection
log.info "Configuration debug info:"
log.info "- Profile: ${workflow.profile}"
log.info "- Config files: ${workflow.configFiles}"
log.info "- Custom config detected: ${hasConfig}"
log.info "- Container engine: ${workflow.containerEngine}"

workflow {
    // Check session configuration details
    log.info "Session config keys: ${session.config.keySet()}"
}
```

### 5.3. Color Support Issues

**Problem**: Colors not displaying correctly in terminal

**Solutions**:

```nextflow title="color_troubleshooting.nf"
include { checkProfileProvided } from 'plugin/nf-core-utils'

// Test color support detection
def supportsColors = System.getenv('TERM') != null &&
                    System.getenv('TERM') != 'dumb' &&
                    System.console() != null

log.info "Color support detected: ${supportsColors}"

// Use appropriate color settings
checkProfileProvided(args, !supportsColors)
```

## 6. Best Practices Summary

### 6.1. Validation Checklist

When implementing configuration validation in your pipeline:

- [ ] Call validation functions early, before heavy processing
- [ ] Provide helpful error messages with solutions
- [ ] Enable colors for better user experience when supported
- [ ] Handle validation failures gracefully
- [ ] Test with various configuration scenarios

### 6.2. Integration Pattern

Follow this standard pattern for nf-core pipelines:

```nextflow title="recommended_pattern.nf" hl_lines="5-7"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// 1. Validate configuration first
include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

// 2. Set up parameters and other imports
params.input = null
params.outdir = "./results"

// 3. Continue with pipeline logic
workflow {
    if (!params.input) {
        error "Please provide --input"
    }

    log.info "Configuration validation passed - starting analysis"
    // Your pipeline processes here...
}
```

### 6.3. Error Prevention Strategy

| User Scenario        | Validation Response               | User Experience               |
| -------------------- | --------------------------------- | ----------------------------- |
| No custom config     | Helpful warning with suggestions  | Guided toward best practices  |
| Profile syntax error | Clear error message with examples | Quick problem resolution      |
| Positional arguments | Warning about ignored parameters  | Prevention of silent failures |
| Good configuration   | Silent validation                 | Smooth pipeline execution     |

## 7. Takeaway

Configuration validation is essential for robust nf-core pipelines:

1. **`checkConfigProvided()`** encourages custom configuration usage for optimal performance
2. **`checkProfileProvided()`** prevents common profile syntax errors and provides helpful guidance
3. **Early validation** catches problems before resource-intensive processing begins
4. **Color-coded messages** improve user experience when supported

The utility transforms potential user frustration into guided learning experiences, making your pipelines more user-friendly and reliable.

## 8. What's Next?

- Explore **[NfcoreNotificationUtils](NfcoreNotificationUtils.md)** for completion notifications and user feedback
- Learn about **[NfCore Utilities](../NfCoreUtilities.md)** for comprehensive pipeline utilities
- Check out **[NextflowPipelineExtension](../NextflowPipelineExtension.md)** for core pipeline functions

:::tip{title="Advanced Usage"}
For complex pipelines with multiple environments, consider implementing environment-specific validation logic using the patterns shown in section 4.2.
:::
