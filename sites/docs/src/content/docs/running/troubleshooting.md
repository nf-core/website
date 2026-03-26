---
title: Troubleshooting
subtitle: Troubleshoot common pipeline issues
shortTitle: Troubleshooting
---

This page covers common issues you may encounter when running nf-core pipelines.

## Troubleshooting basics

### Examine log files

When a pipeline fails, always check `.nextflow.log` for error messages. The work directory for each process contains additional diagnostic files:

- `command.sh`: The command that was executed
- `command.out`: Standard output from the tool
- `command.err`: Standard error from the tool
- `command.log`: Combined output and error
- `exitcode`: The exit code when the job ended
- `command.trace`: Resource usage logs

### Categorize the error

The timing of a failure helps narrow down its cause:

- **Before the first process**: Often caused by an outdated version of Nextflow. An example is the error `Unknown config attribute: projectDir`.
- **During the first process**: Typically a software dependency issue or missing command in `$PATH`.
- **During pipeline execution**: Check the process log files in the work directory.
- **While generating outputs**: Look for a `Missing output file` message followed by tool-specific errors.

## Input and output

### Sample sheet path errors

If the pipeline fails to find input files, verify that all file paths in your sample sheet are correct. A common mistake is trailing spaces at the end of file paths.

### Local iGenomes path overridden

If you are using a local copy of iGenomes but the pipeline defaults to the S3 locations, specify the path directly on the command line:

```bash
nextflow run nf-core/<pipeline> --igenomes_base /path/to/igenomes
```

The `params.igenomes_base` setting in custom config files can be overridden by the load order of configuration files. Using the command-line flag or a `-params-file` ensures the correct path is used.

### Glob pattern not matching files

If only one file is processed when you expect multiple, the shell may be expanding your glob pattern before passing it to Nextflow. Enclose glob patterns in single quotes to prevent shell expansion:

```bash
nextflow run nf-core/<pipeline> --input 'my_*_{R1,R2}.fastq.gz'
```

## Pipeline crashes

### Tool not found at startup

If the pipeline fails immediately with a `command not found` error, you likely have not specified a container or environment profile. Without `-profile`, Nextflow searches for tools in your system `$PATH`. Add a profile to your command:

```bash
nextflow run nf-core/<pipeline> -profile docker
```

Use `docker`, `singularity`, `podman`, or `conda` depending on what is available on your system.

### Stuck on revision

If Nextflow warns that `Project nf-core/<pipeline> is currently stuck on revision`, multiple versions are installed locally and Nextflow cannot determine which one to run. Specify the version explicitly with `-r`:

```bash
nextflow run nf-core/<pipeline> -r 2.1.0 --input ...
```

### Out of memory errors

Exit codes 104, 134, 137, 139, 143, and 247 typically indicate out-of-memory errors. nf-core pipelines automatically retry failed processes with increased resources (default: 3 retries). If retries are exhausted, create a configuration file that adjusts the resource limits for your system and resume:

```bash
nextflow run nf-core/<pipeline> -c custom.config -resume
```

### Steps not executed

If an expected step is skipped, check for typos in optional parameter names. If an upstream process was disabled, it may not have produced the input files required by the skipped step.

## Pipeline updates

### Update not applying

If pulling the latest version of a pipeline does not apply expected fixes, Nextflow may believe you already have the corrected version. Clear the local cache and re-pull:

```bash
rm -r ~/.nextflow/assets/nf-core/<pipeline>
nextflow pull nf-core/<pipeline> -latest
```

## Session errors

### Unable to acquire lock

If you see `Unable to acquire lock on session with ID ...`, a previous Nextflow run in the same directory did not terminate cleanly. Remove the `work` directory to clear the lock:

```bash
rm -rf work
```

To avoid this issue, use `ctrl + c` (not `ctrl + z`) to stop a Nextflow run. `ctrl + z` suspends the process without allowing Nextflow to shut down cleanly.

## Docker errors

### Permission errors when writing to home directory

If a specific process fails with permission errors while writing to `$HOME`, the issue is caused by the `-u $(id -u):$(id -g)` option in the Docker profile. This option emulates your user inside the container but provides no home directory. Rather than disabling this globally, override `containerOptions` for the affected process in your `modules.config`:

```groovy
process {
  withName: TOOL_NAME {
    containerOptions = ''
  }
}
```

Note: `containerOptions` is not supported by Kubernetes or Google Life Sciences executors.

## Network errors

### IPv6 connectivity issues

If you see `WARN: Cannot read project manifest -- Cause: Network is unreachable (connect failed)` on a machine using IPv6, Nextflow and Java default to IPv4. Configure Java to prefer IPv6 by setting this environment variable before running:

```bash
export NXF_OPTS="-Djava.net.preferIPv4Stack=false -Djava.net.preferIPv6Addresses=true"
```

Add this line to `~/.bashrc` or `~/.zshrc` for a permanent fix.

### Docker network isolation on IPv6 machines

If processes inside Docker containers fail to fetch data on IPv6 machines, enable host network mode for the affected process in a custom config file:

```groovy
process {
  withName: YOUR_MODULE {
    containerOptions = '--network host'
  }
}
```

Then run with:

```bash
nextflow run nf-core/<pipeline> -profile docker -c custom.config
```
