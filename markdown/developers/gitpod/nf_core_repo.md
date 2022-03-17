
## Trying your first Gitpod environment

You can run Gitpod with the any nf-core pipeline repository. 

For example, for nf-core RNA-seq pipeline, simply click the green gitpod button or add the Gitpod prefix before the git URL (of https://gitpod.io/# ; to become: https://gitpod.io/#https://github.com/nf-core/rnaseq, 

Once Gitpod has loaded, including the container with all the tools we need, we can go to the terminal and type the following to start the nf-core workflow:

```console
	nextflow run nf-core/rnaseq \
	-profile test,docker \
    --outdir my_result
```

This should run the test data through nf-core rnaseq, using docker with your results in the folder: "my_result". This will take some time to complete.

Using this Gitpod method makes it easy to run and test nf-core pipelines quickly, but lacks the parallelization required to run real datasets.
