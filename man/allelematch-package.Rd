\name{allelematch-package}
\alias{allelematch-package}
\alias{allelematch-tutorial}
\alias{allelematch}
\docType{package}

\title{
	Identification of unique multilocus genotypes
}

\description{
	Tools for the identification of unique multilocus genotypes when both genotyping error and
	missing data may be present; targeted for use with large datasets and databases containing
	multiple samples of each individual (a common situation in conservation genetics, particularly
	in non-invasive wildlife sampling applications). Functions explicitly incorporate missing data
	and can tolerate allele mismatches created by genotyping error.
}

\details{
	\tabular{ll}{
		Package: \tab allelematch\cr
		Type: \tab Package\cr
		Version: \tab 3.0.0\cr
		Date: \tab 2026-07-23\cr
		License: \tab GPL (>= 2)\cr
		Requires: \tab dynamicTreeCut\cr
		LazyLoad: \tab yes\cr
	}

	\subsection{Core Functions}{
	\itemize{
		\item \code{\link{amDataset}}: Produces an input dataset object for allelematch routines.
		\item \code{\link{amPreCheck}}: Pre-screens an amDataset object for severe missing data loads to prevent runtime sorting failures.
		\item \code{\link{amMatrix}}: Produces a dissimilarity matrix.
		\item \code{\link{amPairwise}}: Performs pairwise matching of genotypes.
		\item \code{\link{amCluster}}: Clusters genotypes using dynamic tree cutting.
		\item \code{\link{amAlleleFreq}}: Produces allele frequencies from an amDataset object.
		\item \code{\link{amUnique}}: Identifies unique genotypes.
		\item \code{\link{amUniqueProfile}}: Finds optimal parameters for amUnique.
	}
}

The complete user manual and analytical tutorials describing the operation of the software
in detail are included locally as a package vignette. It is installed automatically with
the package and can be accessed within an active R session using \code{vignette("allelematch")}.\cr

Simulations examining the performance of these tools have also been performed, and results are
available in the publication associated with this package. Please refer to the publication:\cr

Galpern P, Manseau, M, Hettinga P, Smith K, and Wilson P. (2012) allelematch: an R package for
identifying unique multilocus genotypes where genotype error and missing data may be present.
\emph{Molecular Ecology Resources} 12:771-778.\cr

Use \code{citation("allelematch")} for the full citation. Please also use this publication when
citing the package.\cr

An important core element of the package is dynamic tree cutting, and this is made possible via the
cutreeHybrid function within the dynamicTreeCut package for R (Langfelder et al., 2008).
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com}), Todd Cross (\email{todd.cross@gmail.com})
}
