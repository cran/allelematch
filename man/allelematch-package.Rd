\name{allelematch-package}
\alias{allelematch-package}
\alias{allelematch-tutorial}
\alias{allelematch}
\docType{package}
\title{
Identification of unique multilocus genotypes
}
\description{
This package provides tools for the identification of unique of multilocus genotypes when both genotyping error and missing data are present. 
The package is targeted at those working with large datasets and databases containing multiple samples of each individual,
a situation that is common in conservation genetics, and particularly in non-invasive wildlife sampling applications. Functions explicitly incorporate missing data, and can tolerate allele mismatches created by genotyping error.
}
\details{
\tabular{ll}{
Package: \tab allelematch\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2011-09-12\cr
License: \tab GPL-2\cr
Requires: \tab dynamicTreeCut\cr
LazyLoad: \tab yes\cr
}

Supplementary documentation illustrating the use of the software is available as a vignette.  It is installed with the package and linked from the package index help page. An online version is also available at \url{http://nricaribou.cc.umanitoba.ca/allelematch/}.\cr\cr  
Simulations examining the performance of these tools have also been performed, and results are available in the publication associated with this package.\cr\cr
An important core element of the package is dynamic tree cutting, and this is made possible using the
dynamicTreeCut package for R (Langfelder et al., 2008). 
}
\author{
Paul Galpern (\email{pgalpern@gmail.com})
}
\references{
Langfelder P, Zhang B, Horvath S (2008) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. \emph{Bioinformatics}, 24, 719
}
