---
title: Final family representative sequences are faulty in proteinfamilies pipeline 1.3.0
subtitle: Bug inside the extract_family_reps.py script of v1.3.0
category: [pipelines]
type: known_regression
severity: high
publishedDate: "2025-09-29"
reporter:
  - vagkaratzas
reviewer:
  - ewels
  - mashehu
pipelines:
  - name: proteinfamilies
    versions: ["1.3.0"]
modules:
  - EXTRACT_FAMILY_REPS
subworkflows:
  - UPDATE_FAMILIES
  - REMOVE_REDUNDANCY
configuration:
nextflowVersions:
nextflowExecutors:
softwareDependencies:
references:
---

# Issue

The python script `extract_family_reps.py` in this version contained a bug where all sequences from a fasta file would be concatenated to become the family representative sequence, instead of just getting the first sequence of the file.

# Resolution

The relevant python function was updated to correctly extract just the first sequence per family.
The fixed version is available in version 1.3.1 onwards.
