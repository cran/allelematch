\name{amMatrix}
\alias{amMatrix}

\title{
	Produce a dissimilarity matrix for pairs of multilocus genotypes
}

\description{
	Given an amDataset object find the dissimilarities between pairs of multilocus genotypes, taking
	missing data into account.
}

\usage{
	amMatrix(
		amDatasetFocal, 
		missingMethod = 2
		)
}

\arguments{
	\item{amDatasetFocal}{
		An \code{amDataset} object. See \code{\link{amDataset}}.
	}
	\item{missingMethod}{
		The method used to determine the similarity of multilocus genotypes when data is missing.\cr
		The default, (\code{missingMethod = 2}), is preferable in all cases.\cr
		See \code{\link{amMatrix}}.
	}
}

\details{
	This function is the behind-the-scenes workhorse of AlleleMatch, and typically will not be
	called by the user.\cr
	
	\code{missingMethod=2} is the recommended value, and the default, as it has performed
 better in simulations. In this method, missing data matches perfectly with missing data,
 while missing data matches partially with non-missing data.\cr
 
	\code{missingMethod = 1} is retained for experimental purposes. Here, missing data
 matches partially with missing and non-missing data. 
}

\value{
	A distance/dissimilarity matrix of S3 class \code{amMatrix}.
}

\references{
	For a complete vignette, please access via the Data S1 Supplementary documentation and tutorials
	(PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amPairwise}}, \code{\link{amUnique}}
}

\examples{
	\dontrun{
	data("amExample1")
	
	## Produce amDataset object
	myDataset <- 
		amDataset(
			amExample1, 
			missingCode = "-99", 
			indexColumn = 1, 
			metaDataColumn = 2
			)
	
	## Produce dissimilarity matrix
	dissimMatrix <- amMatrix(myDataset)
	}
}