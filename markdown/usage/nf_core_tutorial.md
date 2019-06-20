# _nf-core_ Tutorial

> Originally written for the Nextflow Camp 2019, Barcelona 2019-09-19: ***"Getting started with nf-core"*** _(see [programme](https://www.nextflow.io/nfcamp/2019/phil2.html))_
>
> Duration: **1hr 50**

> Last updated: June 2019

## Table of Contents
* Exercises:
    [1](#exercise-1),
    [2](#exercise-2),
    [3](#exercise-3),
    [4](#exercise-4),
    [5](#exercise-5),
    [6](#exercise-6)
* [Abstract](#abstract)
* [Introduction](#introduction)
* [Installing the _nf-core_ helper tools](#installing-the-nf-core-helper-tools)
* [Listing available _nf-core_ pipelines](#listing-available-nf-core-pipelines)
* [Running _nf-core_ pipelines](#running-nf-core-pipelines)
* [Creating _nf-core_ pipelines](#creating-nf-core-pipelines)
* [Testing _nf-core_ pipelines](#testing-nf-core-pipelines)
* [Releasing _nf-core_ pipelines](#releasing-nf-core-pipelines)

## Abstract

The _nf-core_ community provides a range of tools to help new users get to grips with nextflow - both by providing complete pipelines that can be used out of the box, and also by helping developers with best practices. Companion tools can create a bare-bones pipeline from a template scattered with `TO-DO` pointers and CI with linting tools check code quality. Guidelines and documentation help to get nextflow newbies on their feet in no time. Best of all, the _nf-core_ community is always on hand to help.

In this tutorial we discuss the best-practice guidelines developed by the _nf-core_ community, why they're important and give insight into the best tips and tricks for budding nextflow pipeline developers. ✨

## Introduction
### What is _nf-core_
### What this tutorial will cover
### Where to get help
## Installing the _nf-core_ helper tools
###  Exercise 1
* Install nf-core/tools
* Use the help flag to list the available commands
## Listing available _nf-core_ pipelines
###  Exercise 2
* Use the help flag to print the list command usage
* List all pipelines
* Sort pipelines alphabetically, then by popularity
* Fetch one of the pipelines using `nextflow pull`
* Use `nf-core list` to see if the pipeline you pulled is up to date
* Filter pipelines for those that work with RNA
* Save these pipeline details to a JSON file
## Running _nf-core_ pipelines
### Software requirements for _nf-core_ pipelines
### Usage instructions and documentation
### Running pipelines with test data
### Using _nf-core_ pipelines offline
### Exercise 3
* Install other dependencies (nextflow, docker)
* Print the command-line usage instructions for the nf-core/rnaseq pipeline
* In a new directory, run the nf-core/rnaseq pipeline with the provided test data
* Download the nf-core/rnaseq pipeline for offline use using the `nf-core download` command
## Creating _nf-core_ pipelines
### Using the nf-core template
### `TODO` statements
### How _nf-core_ software packaging works
### Exercise 4
* Make a new pipeline using the template
* Update the readme file to fill in the `TODO` statements
* Add a new process to the pipeline in `main.nf`
* Add the new software dependencies from this process in to `environment.yaml`
## Testing _nf-core_ pipelines
###  Linting _nf-core_ pipelines
### Choosing test data
### nf-core/test_datasets
### Travis CI
### Exercise 5
* Run `nf-core lint` on your pipeline and make note of any test warnings / failures
* Read up on one or two of the linting rules on the nf-core website and see if you can fix some.
* Take a look at `conf/test.config` and switch the test data for another dataset on nf-core/test_data
## Releasing _nf-core_ pipelines
### Forking to _nf-core_
### Initial community review
### Making the first release
## Updating _nf-core_ pipelines
### Pull-requests and branching
### Template updates
### Pipeline version numbers
### Review releases
### Exercise 6
* Use `nf-core bump-version` to update the required version of Nextflow in your pipeline
* Bump your pipeline's version to 1.0, ready for its first release!
* Make sure that you're signed up to the _nf-core_ slack (get an invite on [nf-co.re](https://nf-cor.re)) and drop us a line about your latest and greatest pipeline plans!
* Ask to be a member of the _nf-core_ GitHub organisation by commenting on [this GitHub issue](https://github.com/nf-core/nf-co.re/issues/3)
* If you're a twitter user, make sure to follow the [@nf_core](https://twitter.com/nf_core) account
