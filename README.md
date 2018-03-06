# branchpointer: Prediction of intronic splicing branchpoints

Citation:
Bethany Signal, Brian S Gloss, Marcel E Dinger, Tim R Mercer; Machine learning annotation of human branchpoints, Bioinformatics, btx688, https://doi.org/10.1093/bioinformatics/btx688

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

**Version**: 1.4.1

**Date**: 2018-03-06

**Author**: Beth Signal

**Maintainer**: Beth Signal <b.signal@garvan.org.au>

**Description**: Predicts branchpoint probability for sites in intronic branchpoint windows. 
queries can be supplied as intronic regions; or to evaluate the effects of mutations, SNPs.
