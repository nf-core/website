
## Adding your institution to the contributor's list


In case you couldn't find your organization / group of contributors on the [community pages](https://nf-co.re/community), please add yourself! Just a few lines of YAML and some icon / logos needed.



1. Fork the repository [nf-core/nf-core](https://github.com/nf-core/nf-co.re) to your own github account. How to fork a GitHub project can be found [here](https://guides.github.com/activities/forking/).

2. In your own fork, modify the YAML file nf-core-contributors.yaml in the editor of your choice.

3. Fill in the details about you and your institution directly into the YAML file.
    * Full Name of the Institution
    * Short Name of the Institution (e.g. an acronym)
    * Description (see YAML for examples)
    * Post Address
    * URL website of the Institution
    * Affiliation: could be the subunit of your institution where you are belonging to (if Fancy University is the Institution, then Affilation would be the McFancyDepartment)
    * Affiliation URL
    * Image_fn: filename with the logo (see the point 4.)
    * Contact: Your Name
    * Contact_email: Your E-Mail
    * Contact_github: Your GitHub Shortname
    * Location in Geocoordinates: [longitute, latitude]
    * Twitter: URL of the institutional twitter account if available

4. Furthermore, submit two versions of the institutional logo in the SVG format which you have described it in the YAML file under image_fn:
    1. upload a "white" version to the folder public_html/assets/img/contributors-white. Must be a single monochrome shape with no background colour,
    2. add a "colour" version to the folder public_html/assets/img/contributors-colour.

    * The difference between the white and colour version of the images can be seen by looking at the examples: [white] / [colour] (!!! missing links).
    * The both images should have the same name. E.g., if you name image_fn: foobar.svg, then the both files should be named foobar.svg in the respective folders.
    * If you have only raster images available, please search for a SVG version.
    (_A tip_: wikipedia often uses the SVG format of displayed logos available for the download.)
    * If you were not successfull in this point, skip it and let us know about it in the pull request.
    * _Note:_ Please make sure that you have the permission using logos of your institution in an open source project. The organizations do not usually have any concerns, just want to be notified about it.

5. After you did the hardest part, please open a pull request from your fork. The modifications shall be compared to the master branch [nf-core/nf-core](https://github.com/nf-core/nf-co.re).
The changes have to be reviewed by a member of our core team. Do not hesitate to open the PR. We will provide you the feedback if anything looks weird or hasn't been properly done.



