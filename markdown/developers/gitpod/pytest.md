## How to test a module with pytest in Gitpod

If you enter the gitpod environment for modules (https://gitpod.io/#https://github.com/nf-core/modules), you can run the pytest function in order to debug a pipeline.

You can learn more about `pytest` [here](https://nf-co.re/events/2021/bytesize-17-pytest-workflow)

Once you are in the environment (by clicking the previous link), try runnning an example pytest for an exisiting module:

```console
PROFILE=docker pytest --tag <module_name> --symlink --keep-workflow-wd --git-aware
```