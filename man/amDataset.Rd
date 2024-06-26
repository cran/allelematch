\name{amDataset}
\alias{amDataset}
\alias{print.amDataset}

\title{
	Prepare a dataset for use with allelematch
}

\description{
	Given an input matrix or \code{data.frame} produce a amDataset object suitable for use with
	other allelematch functions.
}

\usage{
	amDataset(
		multilocusDataset, 
		missingCode = "-99", 
		indexColumn = NULL, 
		metaDataColumn = NULL, 
		ignoreColumn = NULL
		)

	\method{print}{amDataset}(x, ...)
}

\arguments{
	\item{multilocusDataset}{
		A \code{matrix} or \code{data.frame} containing samples in rows and alleles in columns.\cr
		Sampling IDs and meta-data may be specified in up to two additional columns.
	}
	\item{missingCode}{
		A character string giving the code used for missing data.\cr
		Missing data may also be represented as \code{NA}.
	}
	\item{indexColumn}{
		Optional.\cr
		A character string giving the column name, or an integer giving the column number containing
		the sampling ID or index information.\cr
		If an index is not supplied, then the function creates an alphabetical index.
	}
	\item{metaDataColumn}{
		Optional.\cr
		A character string giving the column name, or an integer giving the column number containing
		the meta-data.
	}
	\item{ignoreColumn}{
		Optional.\cr
		A vector of character string(s) giving the column name(s) or integer(s) giving the column
		number(s) that should be removed from the input dataset (i.e., those that matching and
		clustering should not consider).
	}
	\item{x}{
		An amDataset object.
	}
	\item{...}{
		Additional arguments to summary.
	}
}

\details{
	Examine \code{\link{amExampleData}} for an example of a typical input dataset in the diploid
	case. (Typically these files will be the CSV output from allele calling software). Sample index
	or ID information and sample meta-data may be specified in two additional columns. Columns can
	optionally be given names, and these are carried through analyses. If column names are not
	given, appropriate names are produced.\cr
 
	Each datum is treated as a character string in \code{allelematch} functions, enabling the mixing
	of numeric and alphanumeric data.\cr

	The multilocus dataset can contain any number of diploid or haploid markers, and these can be in
	any order. Thus in the diploid case there should be two columns for each locus (named, say,
	locus1a and locus1b). Please note that \code{AlleleMatch} functions pay no attention to
	genetics. In other words, each column is considered a comparable state. Thus matching and
	clustering of multilocus genotypes is done on the basis of superficial similarity of the data
	matrix rows, rather than on any appreciation of the allelic states at each locus. See
	\code{\link{amPairwise}} for more discussion.\cr

	For this reason it is important when working with diploid data to ensure that identical
	individuals will have identical alleles in each column. This can be achieved by sorting each
	locus so that in each case the lower length allele appears in, say, a column "locus1a" and the
	higher in column "locus1b." This pattern is likely the default in allele calling software and
	sorting will typically not be required unless data are derived from an unusual source.\cr

	Only one meta-data column is possible with \code{allelematch}. If multiple columns must be
	associated with a given sample for downstream analyses, try pasting them together into one
	string with an appropriate separator, and separating them later when allelematch analyses are
	concluded.
}

\value{
	An \code{amDataset} object.
}

\references{
	For a complete vignette, please access via the Data S1 Supplementary documentation and tutorials
	(PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amPairwise}}, \code{\link{amUnique}}, \code{\link{amExampleData}}
}

\examples{
	\dontrun{
	data("amExample5")
	
	## Typical usage
	myDataset <- 
		amDataset(
			amExample5, 
			missingCode = "-99", 
			indexColumn = 1, 
			metaDataColumn = 2, 
			ignoreColumn = "gender"
			)
	
	## Access elements of amDataset object
	myMetaData <- myDataset$metaData
	mySamplingID <- myDataset$index
	myAlleles <- myDataset$multilocus
	
	## View the structure of amDataset object
	unclass(myDataset)
	}
}