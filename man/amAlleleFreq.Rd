\name{amAlleleFreq}
\alias{amAlleleFreq}
\alias{print.amAlleleFreq}

\title{
	Determine allele frequencies
}

\description{
	Determines the allele frequencies for each locus in a multilocus genetic dataset by 
	first removing missing observations. It requires an \code{amDataset} object and a map 
	vector relating each column of the dataset to a genetic locus.
}

\usage{
	amAlleleFreq(
		amDatasetFocal, 
		multilocusMap = NULL
		)
		
	\method{print}{amAlleleFreq}(x, ...)
}

\arguments{
	\item{amDatasetFocal}{
		An \code{amDataset} object
	}
	\item{multilocusMap}{
		Optionally, a vector of integers or strings giving the mappings onto loci for all 
		genotype columns in \code{amDatasetFocal}.\cr
		When omitted, columns are assumed to be paired (i.e., diploid loci with alleles in 
		adjacent columns).\cr
		See details.
	}
	\item{x}{
		An \code{amDataset} object.
	}
	\item{...}{
		Additional arguments to summary
	}
}

\details{
	This function is called by \code{amUnique}.\cr
	
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
}

\value{
    An \code{amAlleleFreq} object
}

\references{
	For more information and for a complete vignette, please access via the Data S1	
	Supplementary documentation and tutorials (PDF) located at
	<doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amUnique}}
}

\examples{
	\dontrun{
	data("amExample5")
	
	## Produce amDataset object
	myDataset1 <- 
		amDataset(
			amExample5, 
			missingCode = "-99", 
			indexColumn = 1, 
			metaDataColumn = 2
			)
	
	## Usage
	myAlleleFreq <- 
		amAlleleFreq(
			myDataset1, 
			multilocusMap = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 
			11, 11)
			)
	
	## Produce amDataset object, but remove gender column
	myDataset2 <- 
		amDataset(
			amExample5, 
			missingCode = "-99", 
			indexColumn = 1,
			metaDataColumn = 2, 
			ignoreColumn = "gender"
			)
	
	## Because all columns are paired, usage is simpler
	myAlleleFreq <- amAlleleFreq(myDataset2)
	}
}

