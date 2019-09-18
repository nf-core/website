# nf-core/nf-co.re - Contributing Guidelines

Hi there! Many thanks for taking an interest in improving the nf-core website.

## Contribution workflow
If you'd like to write some code for nf-core/nf-co.re, the standard workflow
is as follows:

1. Check that there isn't already an issue about your idea in the
   [nf-core/nf-co.re issues](https://github.com/nf-core/nf-co.re/issues) to avoid
   duplicating work.
    * If there isn't one already, please create one so that others know you're working on this
2. Fork the [nf-core/nf-co.re repository](https://github.com/nf-core/nf-co.re) to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request against the `master` branch and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [basic docs from GitHub](https://help.github.com/articles/fork-a-repo/) or even their [excellent interactive tutorial](https://try.github.io/).


## Tests
When you create a pull request with changes, [Travis CI](https://travis-ci.org/) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.
For now, the only test is for Markdown syntax, using the `markdownlint` package.

## Getting help
For further information/help, please consult the [nf-core/nf-co.re documentation](https://github.com/nf-core/nf-co.re#documentation) and don't hesitate to get in touch on the nf-core `tools` channel on [Slack](https://nf-co.re/join/slack/).
