# branchpointer

##Prediction of intronic splicing branchpoints

###Introduction

###Installation

R-package branchpointer can be installed:

    library(devtools)
    install_github("betsig/branchpointer")

After installation, the package can be loaded into R.

    library(branchpointer)

For details of how to use this package, please see the vignette.

###Dependencies

branchpointer requires [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) to be installed

**Package**: branchpointer

**Type**: Package

**Title**: Prediction of intronic splicing branchpoints

**Version**: 0.3

**Date**: 2016-08-19

**Author**: Beth Signal

**Maintainer**: Beth Signal <b.signal@garvan.org.au>

**Description**: Predicts branchpoint probability for sites in intronic branchpoint windows. 
queries can be supplied as intronic regions; or to evaluate the effects of mutations, SNPs.