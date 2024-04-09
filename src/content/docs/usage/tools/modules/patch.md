---
title: Patch a module
subtitle: Create a patch file for a module
weight: 60
---

If you want to make a minor change to a locally installed module but still keep it up date with the remote version, you can create a patch file using `nf-core modules patch`.

<!-- RICH-CODEX
working_dir: tmp/nf-core-nextbigthing
before_command:  sed "s/process_medium/process_low/g" modules/nf-core/fastqc/main.nf > modules/nf-core/fastqc/main.nf.patch && mv modules/nf-core/fastqc/main.nf.patch modules/nf-core/fastqc/main.nf
-->

![`nf-core modules patch fastqc`](/images/tools/nf-core-modules-patch.svg)

The generated patches work with `nf-core modules update`: when you install a new version of the module, the command tries to apply
the patch automatically. The patch application fails if the new version of the module modifies the same lines as the patch. In this case,
the patch new version is installed but the old patch file is preserved.

When linting a patched module, the linting command will check the validity of the patch. When running other lint tests the patch is applied in reverse, and the original files are linted.
