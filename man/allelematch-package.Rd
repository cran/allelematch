\name{allelematch-package}
\alias{allelematch-package}
\alias{allelematch-tutorial}
\alias{allelematch}
\docType{package}
\title{
Identifying Unique Multilocus Genotypes where Genotyping Error and Missing Data may be Present
}
\description{
This package provides tools for the identification of unique of multilocus genotypes when both genotyping error and missing data may be present. 
The package is targeted at those working with large datasets and databases containing multiple samples of each individual,
a situation that is common in conservation genetics, and particularly in non-invasive wildlife sampling applications.
Functions explicitly incorporate missing data, and can tolerate allele mismatches created by genotyping error.
}
\details{
\tabular{ll}{
Package: \tab allelematch\cr
Type: \tab Package\cr
Version: \tab 2.5.1\cr
Date: \tab 2019-04-08\cr
License: \tab GPL-3\cr
Depends: \tab graphics, stats, utils\cr
Imports: \tab dynamicTreeCut\cr
LazyLoad: \tab yes\cr
}

Supplementary documentation describing the operation of the software in detail and illustrating the use of the software using tutorials is available as a vignette.
It is installed with the package and linked from the package index help page.

Simulations examining the performance of these tools have also been performed, and results are available in the publication associated with this package.  Please refer to the publication:\cr

Galpern P, Manseau, M, Hettinga P, Smith K, and Wilson P.  (2012)  allelematch: an R package for identifying unique multilocus genotypes where genotype error and missing data may be present.
\emph{Molecular Ecology Resources} 12:771-778.\cr

Use \code{citation("allelematch")} for the full citation.  Please also use this publication when citing the package.\cr

An important core element of the package is dynamic tree cutting, and this is made possible using the
dynamicTreeCut package for R (Langfelder et al., 2008). 
}
\author{
Paul Galpern (\email{pgalpern@gmail.com})
}
\references{
Langfelder P, Zhang B, Horvath S. (2008) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. \emph{Bioinformatics}, 24, 719
}
