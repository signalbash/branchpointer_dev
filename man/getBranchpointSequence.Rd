% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictBranchpoints.R
\name{getBranchpointSequence}
\alias{getBranchpointSequence}
\title{Get branchpoint sequence features
Gets intronic sequence covering the branchpoint window and extracts predictive features}
\usage{
getBranchpointSequence(query, uniqueId = "test", queryType,
  workingDirectory = ".", genome = NA, bedtoolsLocation = NA,
  BSgenome = NULL, useParallel = FALSE, cores = 1, rmChr = FALSE)
}
\arguments{
\item{query}{branchpointer query GenomicRanges}

\item{uniqueId}{unique string identifier for intermediate .bed and .fa files.}

\item{queryType}{type of branchpointer query. "SNP" or "region" or "indel".}

\item{workingDirectory}{directory where intermediate .bed and .fa are located}

\item{genome}{.fa genome file location}

\item{bedtoolsLocation}{bedtools binary location (which bedtools)}

\item{BSgenome}{BSgenome object}

\item{useParallel}{use parallelisation to speed up code?}

\item{cores}{number of cores to use in parallelisation (default = \code{1})}

\item{rmChr}{remove "chr" before chromosome names before writing bed file.
Required if genome sequence names do not contain "chr"}
}
\value{
GenomicRanges with all features required to predict branchpoint probability scores
}
\description{
Get branchpoint sequence features
Gets intronic sequence covering the branchpoint window and extracts predictive features
}
\author{
Beth Signal
}
\keyword{internal}
