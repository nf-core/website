---
title: Manual synchronisation
subtitle: How to manually synchronise your nf-core pipeline with the latest nf-core/tools template
---

There are rare cases, when the synchronisation needs to be triggered manually,
i.e. it was not executed during an `nf-core/tools` release on Github, or when you want to perform a targeted sync.

Note that automated PR system is only applicable to official nf-core pipelines, homemade pipelines based on nf-core standards/modules created with `nf-core create` have to be updated following this manual synchronisation procedure.

You can do so by running the `nf-core sync` command:

```bash
cd my_pipeline
git checkout dev # or your most up to date branch
nf-core sync -d .
```

Note that the `sync` command assumes that you have a branch called `TEMPLATE`, so you may need to pull this from the upstream nf-core repository if you are working on a fork:

```bash
git remote add upstream https://github.com/nf-core/PIPELINE.git
git checkout --track upstream/TEMPLATE
```

Remember to go back to your `dev` branch as above before running `nf-core sync`.

Much of the merging process should then be the same as described above with the automated pull requests.

In case of manual synchronisation of a homemade pipeline and if you want to have a PR opened to your `dev` branch, you can use this `nf-core sync` template command:

```bash
nf-core sync \
   --dir [pipeline_dir]
   --from-branch dev \
   --pull-request \
   --username [GitHub_username] \
   --github-repository [GitHub_pipeline_URL]
```
