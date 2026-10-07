---
title: Secrets in the debug log artifact of the AWS test workflows
subtitle: The AWS test workflows in the pipeline template uploaded a debug log containing secrets as a public workflow artifact, letting anyone signed in to GitHub read them
category:
  - pipelines
type:
  - security
severity: high
publishedDate: "2026-10-06"
reporter:
  - ewels
reviewer:
  - mashehu
pipelines:
  - all with template version 2.9 - 2.13.1
  - all with template version 3.4.0 - 4.1.0
modules:
subworkflows:
configuration:
nextflowVersions:
nextflowExecutors:
softwareDependencies:
references:
  - title: seqeralabs/action-tower-launch v2.4.0 release
    description: Always removes the access token from the debug log, including when a launch fails
    url: https://github.com/seqeralabs/action-tower-launch/releases/tag/v2.4.0
  - title: seqeralabs/action-tower-launch v2.3.0 release
    description: Stops logging full HTTP request and response bodies by default
    url: https://github.com/seqeralabs/action-tower-launch/releases/tag/v2.3.0
  - title: Upstream fix for the access token on the failure path
    description: Explains why the token was left in the log when a launch failed
    url: https://github.com/seqeralabs/action-tower-launch/pull/47
---

# Issue

:::note{title="tl;dr"}
The AWS test workflows uploaded a debug log that could contain the Slack credentials used for test notifications. In nf-core/rnaseq, it could also contain the Seqera Platform access token.
:::

The nf-core pipeline template ships two GitHub Actions workflows (`awstest.yml` and `awsfulltest.yml`) that launch test runs on AWS through Seqera Platform with the [`seqeralabs/action-tower-launch`](https://github.com/seqeralabs/action-tower-launch) action.
The action writes a debug log, `tower_action_*.log`, and the template uploads it as a workflow artifact called "Seqera Platform debug log file".
On public repositories, any signed-in GitHub user can download these artifacts for as long as they are kept (90 days by default).

Versions of the action before v2.3.0 ran the Seqera Platform CLI in verbose mode, which writes full HTTP requests into that log.
The launch request contains the pipeline parameters and Nextflow configuration that the workflow passes to the action, including any secrets written into them.
The template's `awsfulltest.yml` passes a Slack credential this way, so that test results get posted to Slack:

- nf-core/tools 2.9 through 2.13.1 and 3.4.0 through 3.5.2: the `MEGATESTS_ALERTS_SLACK_HOOK_URL` webhook URL, in `parameters`
- nf-core/tools 4.0.0 through 4.1.0: the `NFSLACK_BOT_TOKEN` Slack bot token, in `nextflow_config`

nf-core/tools 2.14.0 through 3.3.2 are not affected, because their template looked for the log under the wrong name (`seqera_platform_action_*.log`) and never uploaded it (lucky us 😬).

The action did remove the Seqera Platform access token (`TOWER_ACCESS_TOKEN`) from the log, but only after a successful launch.
When the launch failed, the action stopped before that step and left the token in the file.
The template only uploads the log after a successful launch, so this only matters for pipelines that changed the upload step to also run on failure.
Among nf-core pipelines, only nf-core/rnaseq did this: since April 2026, its `awsfulltest.yml` uploads the log with `if: always()`, so the debug logs of its failed full-size test launches contain the access token.

In short: any signed-in GitHub user could download the debug logs of an affected pipeline and read the Slack credentials, and for nf-core/rnaseq also the Seqera Platform access token.

Pipelines using the `@v2` tag of the action (nf-core/tools 2.9 through 3.5.2) received the fixed version automatically when the tag moved.
Pipelines created or last synchronized with nf-core/tools 4.0.0 through 4.1.0 pin the action to v2.2.0 and are still affected.
The fix touches the following line in `.github/workflows/awstest.yml` and `.github/workflows/awsfulltest.yml`:

```diff title=".github/workflows/awstest|awsfulltest.yml"
-        uses: seqeralabs/action-tower-launch@51565b514bff1827cf34620de25d0055759f1fc9 # v2
+        uses: seqeralabs/action-tower-launch@v2
```

We are going back to the non-hashed version of this actions, since it is maintained by nf-core members and this allows us to push fixes to the action directly.

# Resolution

:::note{title="tl;dr"}
An automated patch PR against the default branch will move nf-core pipelines to the fixed version of the action. It will be merged by the nf-core infrastructure team.
:::

[v2.4.0](https://github.com/seqeralabs/action-tower-launch/releases/tag/v2.4.0) of the action fixes this: it no longer logs request bodies by default, and it always removes the access token from the debug log, even when a launch fails.
The `@v2` tag of the action points to that release.

Pipeline maintainers can merge the automated `patch` pull request opened against their default branch if the infrastructure team isn't fast enough for them.
Releases and the full-size tests on release pull requests both use the workflow file from the default branch, so the fix takes effect for them as soon as it is merged.
The `dev` branch gets the same change with the next template sync.
No new pipeline release is required after the fix is merged, because it only affects CI workflows and not the actual pipeline code.

Pipelines outside nf-core that use `seqeralabs/action-tower-launch` should update it to v2.4.0 or later.
If you pass secrets in `parameters` or `nextflow_config`, you should also check old "Seqera Platform debug log file" artifacts, delete those that contain secrets, and rotate the affected credentials.
