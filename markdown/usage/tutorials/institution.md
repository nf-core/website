
## Adding Your Institution to the Contributor's List
### How To Add Yourself to the Community


In case you couldn't find your organization / group of contributors on the [community pages](https://nf-co.re/community), please add yourself. Just a few lines of YAML and your logos needed.



1. Fork the repository [nf-core/nf-core](https://github.com/nf-core/nf-co.re) to your own GitHub account. Look [here](https://guides.github.com/activities/forking/) for advice how to fork a GitHub project.

2. In your own fork, modify the YAML file 'nf-core-contributors.yaml' in the editor of your choice.

3. Fill in the details about you and your institution directly into the YAML file.
    * Full Name of the Institution
    * Short Name of the Institution (e.g. an acronym)
    * Description (see YAML for the examples)
    * Postal Address
    * Institution URL
    * Affiliation: could be the subunit of your institution where you are belonging to (if Fancy University is the Institution, then Affilation would be the McFancyDepartment)
    * Affiliation URL
    * Image_fn: filename with the logo (see the point [4](#4).)
    * contact: Your Name
    * contact_email: Your E-Mail
    * contact_github: Your GitHub Shortname
    * location in geocoordinates: [longitude, latitude]
    * twitter: URL of the institutional twitter account, if available

4. Furthermore, submit two versions of the institutional logo in SVG format which you have described it in the YAML file under image_fn:
    1. Upload a "white" version to the folder public_html/assets/img/contributors-white. Must be a single monochrome shape with no background colour.
    2. Upload a "colour" version to the folder public_html/assets/img/contributors-colour.

    * These examples in [white](https://github.com/nf-core/nf-co.re/tree/master/public_html/assets/img/contributors-white) / [colour](https://github.com/nf-core/nf-co.re/tree/master/public_html/assets/img/contributors-colour) indicate to you the difference between the versions.
    * Both images should have the same name. I.e., if you name them image_fn: foobar.svg, then both files should be named foobar.svg in the respective folders.
    * If you have only raster images available, please search for a SVG version.
    (_A tip_: wikipedia often uses the SVG format of displayed logos available for the download.)
    * If you were not successful at this point, skip it and let us know about it in the pull request.
    * _Note:_ Please make sure that you have the permission to use logos of your institution in an open source project. Organizations often don't have any concerns about it, just want to be notified.

5. After you did the hardest part, please open a pull request from your fork. The modifications shall be compared to the master branch [nf-core/nf-core](https://github.com/nf-core/nf-co.re).
The changes have to be reviewed by a member of our core team. Do not hesitate to open the PR. We will kindly provide you a feedback if anything looks weird or hasn't been properly done.
