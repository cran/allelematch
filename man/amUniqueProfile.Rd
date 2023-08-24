\name{amUniqueProfile}
\alias{amUniqueProfile}

\title{
	Determine optimal parameter values for the identification of unique genotypes
}

\description{
	Function to automatically run \code{amUnique} at a sequence of parameter values to
	determine an optimal setting, and optionally plot the result
}

\usage{
	amUniqueProfile(
		amDatasetFocal, 
		multilocusMap = NULL, 
		alleleMismatch = NULL, 
		matchThreshold = NULL, 
		cutHeight = NULL, 
		guessOptimum = TRUE, 
		doPlot = TRUE, 
		consensusMethod = 1, 
		verbose = TRUE
		)
}

\arguments{
	\item{amDatasetFocal}{
		An \code{amDataset} object containing genotypes in which an unknown number of
		individuals are sampled multiple times.
	}
	\item{multilocusMap}{
		Optionally a vector of integers or strings giving the mappings onto loci for all
		genotype columns in amDatasetFocal.\cr
		When omitted, columns are assumed to be paired (i.e., diploid loci with alleles in
		adjacent columns).\cr
		See details.
	}
	\item{alleleMismatch}{
		A vector giving a sequence, where elements give the maximum number of mismatching
		alleles which will be tolerated when identifying individuals; also known as the
		m-hat parameter.\cr
		If specified, then \code{matchThreshold} and \code{cutHeight} should be omitted;
		all three parameters are related.\cr
		See \code{\link{amUnique}} for details.
	}
	\item{matchThreshold}{
		A vector giving a sequence, where elements give the minimum dissimilarity score
		which constitutes a match when identifying individuals; also known as the s-hat
		parameter.\cr
		If specified, then \code{alleleMismatch} and \code{cutHeight} should be omitted;
		all three parameters are related.\cr
		See \code{\link{amUnique}} for details.
	}
	\item{cutHeight}{
		A vector giving a sequence, where elements give the \code{cutHeight} parameter
		used in dynamic tree cutting by \code{amCluster}; also known as the d-hat
		parameter.\cr
		If specified, then \code{alleleMismatch} and \code{matchThreshold} should be
		omitted; all three parameters are related.\cr
		See details.
	}
	\item{doPlot}{
		If \code{doPlot = TRUE}, then a plot showing summary data from multiple runs of
		\code{amUnique} is produced.
	}
	\item{guessOptimum}{
		If \code{guessOptimum = TRUE}, then will estimate the optimal value of the
		parameter being profiled by examining the profile for the first minimum associated
		with a drop in multiple matches as sensitivity to differences among samples
		decreases.
	}
	\item{consensusMethod}{
		The method (an integer) used to determine the consensus multilocus genotype from a
		cluster of multilocus genotypes.\cr
		See \code{\link{amCluster}} for details.
	}
	\item{verbose}{
		If \code{verbose = TRUE}, then the progress of the profiling will be reported to
		the console.
	}
}

\details{
	Selecting the appropriate value for \code{alleleMismatch}, \code{cutHeight}, or
	\code{matchThreshold} is an important task. Use this function to assist in this
	process. Typically the optimal value of any of these parameters is found where the
	number of multiple matches is minimized (the majority of samples are similar to only
	one unique genotype). Usually there is a minimum when these parameters are set to be
	very sensitive to differences among samples (i.e., \code{alleleMismatch} or
	\code{cutHeight} are 0, \code{matchThreshold} is 1). Simulations suggest that the next
	most sensitive minimum in multiple matches is the optimal value. This minimum will
	often be associated with a drop in multiple matches as sensitivity drops. For more
	discussion of this important step, see the Data S1 Supplementary documentation and
	tutorials (PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.\cr
	
	Using \code{guessOptimum = TRUE} will attempt to estimate the location of this minimum
	and add it to the profile plot. Manual assessment of this estimate using the plot is
	strongly recommended.\cr
	
	If none of \code{alleleMismatch}, \code{cutHeight}, or \code{matchThreshold} is given,
	the function runs a sequence of values for \code{alleleMismatch} as follows:
	\code{seq(from = 0, to = floor(ncol(amDatasetFocal$multilocus) * 0.4), by = 1)}\cr
	
	\code{multilocusMap} is often not required, as \code{amDataset} objects will typically
	consist of paired columns of genotypes, where each pair is a separate locus. In cases
	where this is not the case (e.g., gender is given in only one column), a map vector
	must be specified.\cr
	
	Example: \code{amDataset} consists of gender followed by 4 diploid loci in paired
	columns\cr
	\code{multilocusMap = c(1, 2, 2, 3, 3, 4, 4, 5, 5)}\cr
	or equally\cr
	\code{multilocusMap=c("GENDER", "LOC1", "LOC1", "LOC2", "LOC2", "LOC3", "LOC4",
	"LOC4")}
	
	For more information on selecting \code{consensusMethod} see \code{\link{amCluster}}.
	The default \code{consensusMethod = 1} is typically adequate.
}

\value{
	A \code{data.frame} containing summary data from multiple runs of \code{amUnique}
}

\references{
	For a complete vignette, please access via the Data S1 Supplementary documentation and
	tutorials (PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amUnique}}
}

\examples{
	\dontrun{
	data("amExample2")
	
	## Produce amDataset object
	myDataset <- 
		amDataset(
			amExample2, 
			missingCode = "-99", 
			indexColumn = 1, 
			metaDataColumn = 2
			)

	## Usage (uncomment)
	myUniqueProfile <- 
		amUniqueProfile(myDataset)
	
	## Data set with gender information
	data("amExample5")
	
	## Produce amDataset object
	myDataset2 <- 
		amDataset(
			amExample5, 
			missingCode = "-99", 
			indexColumn = 1, 
			metaDataColumn = 2
			)
	
	## Usage
	myUniqueProfile <- 
		amUniqueProfile(
			myDataset2,
			multilocusMap = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10,
			11, 11))
	}
}