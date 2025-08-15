---
title: "10. FAQ & Debugging"
subtitle: Best practices, common issues, and solutions for nf-test
weight: 100
---

## Debugging Techniques

### Using println for Debugging

Print statements can help debug test issues. They must go within the `then` block, prior to `assertAll`:

```groovy
then {
    def unstable_patterns_auth = [
        '**/mapped_reads_gc-content_distribution.txt',
        '**/genome_gc_content_per_window.png',
        '**/*.{svg,pdf,html}',
        '*.{svg,pdf,html}',
        '**/DamageProfiler.log',
    ]

    println("unstable_patterns_auth: " + unstable_patterns_auth)

    assertAll(
        { assert snapshot( stable_content_authentication, stable_name_authentication*.name ).match("authentication") },
        // ... more assertions
    )
}
```

### Using Diff Tools for Snapshot Comparison

When snapshot tests fail, nf-test automatically shows differences between expected and actual snapshots to help identify mismatches. By default, nf-test uses the `diff` tool with side-by-side comparison (`-y`) and 200-character width (`-W 200`).

#### Customizing Diff Tool Arguments

You can customize the diff tool behavior using the `NFT_DIFF_ARGS` environment variable:

```bash
export NFT_DIFF_ARGS="<your_custom_arguments>"
nf-test test
```

#### Using Alternative Diff Tools

For better visualization, you can change the diff tool entirely using the `NFT_DIFF` environment variable. For example, to use `icdiff` (improved colored diff):

```bash
# Install icdiff first
pip install icdiff

# Set environment variables
export NFT_DIFF="icdiff"
export NFT_DIFF_ARGS="-N --cols 200 -L expected -L observed -t"

# Run tests
nf-test test
```

#### Benefits of Enhanced Diff Tools

These tools help quickly identify:

- **Text differences**: Line-by-line changes with color highlighting
- **Structural changes**: Better visualization of file organization changes
- **Content mismatches**: Clear indication of what changed between snapshots

**Pro tip**: Use enhanced diff tools like `icdiff` during development to more easily spot snapshot differences and understand why tests are failing.

For more details, see the [official nf-test snapshot differences documentation](https://www.nf-test.com/docs/assertions/snapshots/#snapshot-differences).

## Known Issues and Container Considerations

When using nf-test with Docker, Singularity, or Conda, be aware of environment-specific issues that can cause mismatched hashes:

### Tips for Handling Mismatched Hashes

1. **Consistent Environment**: Ensure consistent environments across containers
2. **Identical Base Images**: Use same base images for Docker/Singularity containers
3. **Pin Software Versions**: Explicitly pin software versions and dependencies
4. **Isolate Non-Deterministic Elements**: Identify and isolate non-deterministic elements
5. **Reproducible Conda Environments**: Use `conda list --explicit` for exact environment recreation
6. **Review Container Caching**: Be cautious with caching mechanisms
7. **Consistent Filesystem Paths**: Ensure path consistency within containers
8. **Regular Updates**: Regularly update and test containers

## Best Practices Summary

1. **Start Simple**: Begin with basic success checks and version snapshots
2. **Capture Everything Possible**: Aim to snapshot complete outputs when stable
3. **Handle Instability Gracefully**: Use selective snapshots for unstable content
4. **Use Meaningful Tags**: Name snapshots descriptively for easier debugging
5. **Debug Systematically**: Use println statements to understand output structure
6. **Review Snapshot Changes**: Always review snapshot file changes in code reviews
7. **Test Edge Cases**: Include tests for error conditions and boundary cases

### Additional Reading

- [nf-test Documentation](https://code.askimed.com/nf-test/docs/getting-started/)
- [Updating Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#updating-snapshots)
- [Cleaning Obsolete Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#cleaning-obsolete-snapshots)
- [Constructing Complex Snapshots](https://code.askimed.com/nf-test/docs/assertions/snapshots/#constructing-complex-snapshots)
