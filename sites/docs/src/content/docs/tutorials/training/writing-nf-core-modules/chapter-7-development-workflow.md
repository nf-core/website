---
title: "Chapter 7: Development workflow"
subtitle: "A summary of the development workflow for contributing to the community"
---

In this chapter, we are going to give you an overview into the development workflow for contributing to the nf-core/modules repository.

As you saw, we ran this training within a fork and clone of the community nf-core modules repository that has more than 1400 already existing modules.
However, everything you have learnt so far, can be used for your own custom modules - whether in your own private repository or even a 'local' module within your pipeline.
But as, nf-core is very pro-community and saving everyone's time by not re-inventing the wheel, we encourage you to contribute your modules to the nf-core/modules repository, meaning that everyone can benefit!

### Check

The first thing to do before you start writing a module is to check that someone has not already made (or is making) it!

To check for the status of the module, you do the following steps:

1. Search for the module the nf-core website on the [modules page](https://nf-co.re/modules)
2. Search for [open PRs](https://github.com/nf-core/modules/pulls) on the nf-core/modules repository
3. Search for [open issues](https://github.com/nf-core/modules/issues) on the nf-core/modules repository

Once you've check these, you can do the determine the next step. So, if:

1. The module already exists, then you most likely can skip any writing and install it directly into your pipeline!
2. Someone is already actively working on writing the module as seen by an open PRm then you can offer to help out (or take over if the development of the module has stalled).
3. There is an issue, then you can either assign yourself or if someone is already assigned - ask if they want any help!

But if the module is not being worked on and no issue exists, [create an issue](https://github.com/nf-core/modules/issues) to let the community know that you are working on it!

Once you've assigned yourself, you should fork and clone the nf-core/modules repository and create a new branch for your module, and follow the steps in the previous chapters to write and test the module.

### Write

To recap, you can use the following command to generate the boilerplate template files for your module:

```bash
nf-core modules create <toolname>/<subcommand>
```

And then go through out each file, filling in the necessary information as guided by the TODO comments, the [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules), and this training material.

### Test

Once you've filled in all the necessary bits and commands, you can test your module and generate the snapshot using the following command:

```bash
nf-core modules test <toolname>/<subcommand>
```

Once you've checked the snapshot file looks good, there is one more step before submitting the module to the official nf-core modules repository.

### Lint

Another very important part of nf-core is standardisation and consistency (as much as we can).
Therefore to help ensure you follow the specifications correctly, the nf-core tools package has a linting tool that can check your module for any issues.

:::note
You do not necessarily have to run linting if you're writing just a custom module for your own use, or for your own module repository.
However it is still recommended to pick up on little issues that might occur.
:::

If you run the following command:

```bash
nf-core modules lint <toolname>/<subcommand>
```

You should get a list of any issues that need to be fixed before submitting the module.
Furthermore, in some cases it can also fix some of issues for you, in which case it will inform you with the relevant command.

Make sure the linting command correctly fixes things for you, in a few cases it may add missing sections you still need to fill in.

:::warning
The nf-core linting command cannot check for all aspects of the [nf-core module specification](https://nf-co.re/docs/guidelines/components/modules)!

Therefore, it is always a good idea to check the specification yourself to ensure you have followed it correctly.
:::

### Submit

If all is good, and you want to share your module with the community, you can commit and push your new module to your fork:

```bash
git commit -am "Add new module <toolname>/<subcommand>"
git push
```

and then open a pull request from your fork to the nf-core/modules repository!

At this point you will need to get a review.
The best way to get a community members attention and get a timely review is to ask on the nf-core slack's [`#request-review`](https://nfcore.slack.com/archives/CQY2U5QU9) channel.
If you're not already on the nf-core slack, you can find the join instructions [here](https://nf-co.re/join).

Once you get your review, you may need to make some corrections.
Don't close your PR - just make the changes on your local clone of the repository and push them to the same branch, and request review (either from the same reviewer, or again from the #request-review channel).

Once you get your ✅, you can merge the module in, and it will be available for everyone to use!

And of course, thank you for contributing to helping make the wider Nextflow community's development experience even smoother!

### Summary

In this chapter we have summarised the development workflow for contributing to the nf-core/modules repository.
We have recapped the most important commands while writing the module ,and then how to lint the module and the procedure for submitting the module to the community repository.

In the next chapter, we will describe how to use the modules in your own Nextflow and official nf-core pipelines.
