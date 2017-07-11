# branchpointer: Prediction of intronic splicing branchpoints

See: http://www.biorxiv.org/content/early/2016/12/14/094003

This is the development version of branchpointer.
For the current BioConductor version see:

https://bioconductor.org/packages/release/bioc/html/branchpointer.html or

https://github.com/Bioconductor-mirror/branchpointer


### Installation

R-package branchpointer can be installed:

    library(devtools)
    install_github("betsig/branchpointer_dev")

After installation, the package can be loaded into R.

    library(branchpointer_dev)

For details of how to use this package, please see the vignette.

### Dependencies

branchpointer requires [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) to be installed if using a local .fa file for sequence retrieval.
Alternatively, an R BSgenome object can now be used (see vignette).

**Package**: branchpointer

**Type**: Package

**Title**: Prediction of intronic splicing branchpoints

**Version**: 1.3.0

**Date**: 2017-07-10

**Author**: Beth Signal

**Maintainer**: Beth Signal <b.signal@garvan.org.au>

**Description**: Predicts branchpoint probability for sites in intronic branchpoint windows. 
queries can be supplied as intronic regions; or to evaluate the effects of mutations, SNPs.
