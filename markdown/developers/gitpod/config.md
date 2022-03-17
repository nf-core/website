## Configuration of Gitpod

Within each git repository, the main file that controls the gitpod environment is the `.gitpod.yml` file, that contains the instructions on which environment to build and which tools to install. 

Check out the nf core `.gitpod.yml` file [here](https://github.com/nf-core/nf-co.re/blob/master/.gitpod.yml). You will often see five main sections:

1. **github** - allows configuration of github. e.g. Allows gitpod to create prebuilds for branches.
2. **vscode** - allows vscode extensions within your environment.
3. **ports**  - opens a port to serve traffic to a public URL
4. **tasks**  - this tells gitpod to run particular jobs. Usually you will see the following:

`- init:` sections can be used to install packages as a pre-build, so it doesn't have to run each time you open an environment.

`- command:` sections execute the given lines of code on every workspace startup.

5. **image** - a container image to pull into Gitpod. Many nf-core pipelines use the image `nfcore/gitpod:latest`. This allows the gitpod environment to contain working nextflow and nf-core scripts and other essential tools such as docker. 

For more detailed information about these settings, check out the extensive docs at Gitpod [here](https://www.gitpod.io/docs/config-gitpod-file).
