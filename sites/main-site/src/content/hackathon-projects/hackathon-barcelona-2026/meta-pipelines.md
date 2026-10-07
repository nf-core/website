---
title: Meta-pipelines and pipeline composition
category: tooling
slack: https://nfcore.slack.com/channels/wg-meta-pipelines
location: Barcelona
leaders:
  pinin4fjords:
    name: Jonathan Manning
---

## Goal

Work out how nf-core pipelines should be made composable into meta-pipelines, and which standards need to be in place first.

This is deliberately open ended. Composition is something nf-core should support, but how to do it is not settled, and part of the session is agreeing the path.

## Description

Many analyses chain several nf-core pipelines, for example [nf-core/rnaseq](https://nf-co.re/rnaseq) followed by [nf-core/differentialabundance](https://nf-co.re/differentialabundance). Today that means running them separately and passing files between them by hand, which loses `-resume` across the whole analysis and adds orchestration glue.

Nextflow [pull request #7213](https://github.com/nextflow-io/nextflow/pull/7213) (see the [pipeline composition ADR](https://github.com/nextflow-io/nextflow/blob/master/adr/20260608-pipeline-composition.md)) lets one pipeline be included in another. [pinin4fjords/rnaseq-diffabundance-meta](https://github.com/pinin4fjords/rnaseq-diffabundance-meta) is a proof of concept that composes rnaseq and differentialabundance into one DAG, one work directory and one `-resume`. The meta-pipeline is just `main.nf`, a config shell and a test profile on top of vendored copies of the two pipelines.

The proof of concept shows that composition is feasible, but nf-core is not ready for it yet. It needed changes to the component pipelines ([rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966), [differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758)) that go beyond current conventions, and it depends on unreleased Nextflow and a patched nf-schema:

- **A typed `params {}` block** as the pipeline interface, with inputs that can be a path on the command line or a channel from an including pipeline
- **Params passed in** rather than read from the global `params`, which in a composed run belongs to the meta-pipeline
- **Config shipped as loadable files**, because an included pipeline's `nextflow.config` is not loaded, and with selectors that work whatever alias the pipeline is included under
- **Tool arguments as process inputs** rather than `ext.args` closures that read `params`, which changes how modules look compared with nf-core/modules
- **Outputs returned as channels** and published by the meta-pipeline, with an `output {}` block in each pipeline
- **Locating files relative to the pipeline** (`moduleDir`) rather than `projectDir` or `bin/`

Several of these depend on work that is itself still being agreed: records and typed pipelines (the rnaseq output-records work in [#1945](https://github.com/nf-core/rnaseq/pull/1945)), the static type system, the params block and workflow outputs, and what `nf-core pipelines lint` and the template should do with them ([nf-core/tools#4491](https://github.com/nf-core/tools/issues/4491), [#4493](https://github.com/nf-core/tools/issues/4493)). A sensible order of work may well be to settle those first and build composition on top.

## Tasks

This is a discussion-led session. Possible directions, to be chosen on the day:

- Walk through the proof of concept, and try composing other pairs of pipelines to see what breaks
- List the blockers: which standards (records, static types, the params block, workflow outputs, module input conventions) must exist before composition is worth standardising, and what state each is in
- Propose an order of work, so that each standard lands before the composition changes that depend on it
- Discuss how tool arguments should reach modules, and what that means for `ext.args` and for sharing modules with nf-core/modules
- Discuss what the template and `nf-core pipelines lint` need to support typed, composable pipelines
- Discuss how a meta-pipeline should pin, vendor and update its component pipelines, and how it should be tested
- Agree next steps, and write up the outcome as guidance or an RFC
