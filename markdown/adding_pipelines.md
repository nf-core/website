# How to add a new pipeline

1. Join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3)
2. Create a pipeline repository in the organisation
    * _If starting from scratch_
        * Ask an admin to create a new pipeline repository for you and add you as a collaborator.
        * This is the best option, as this repository will then be recognised as the "head fork" by GitHub.
    * _If you already have a pipeline_
        * Talk to us about it on the [nf-core gitter chat](https://gitter.im/nf-core/Lobby)
        * Fork your pipeline to the nf-core GitHub organisation
        * Go through [this guide](sync#setup) in order to set up the template syncing
3. Make sure that your pipeline `README.md` file has a big warning on the front saying that it's under development
4. Work on the `dev` branch until the automatic lint tests are passing
5. When you're happy, ping [@nf-core/core](https://github.com/orgs/nf-core/teams/core) for a code review
5. Once the pipeline has been approved, you can remove the development warning on the main readme and the pipeline. Tag a release and the website will be updated.
