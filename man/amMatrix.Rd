\name{amMatrix}
\alias{amMatrix}
\title{
Produce a dissimilarity matrix for pairs of multilocus genotypes
}
\description{
Given an amDataset object find the dissimilarities between pairs of multilocus genotypes,
taking missing data into account.
}
\usage{
amMatrix(amDatasetFocal, missingMethod = 2)
}
\arguments{
  \item{amDatasetFocal}{
  An amDataset object.  See \code{\link{amDataset}}.
}
  \item{missingMethod}{
  The method used to determine the similarity of multilocus genotypes when data is
  missing.  The default (=2) is preferable in all cases.  Please see details.
}
}
\details{
 This function is the behind-the-scenes workhorse of AlleleMatch, and typically
 will not be called by the user.\cr
 
 \code{missingMethod=2} is the recommended value, and the default, as it has performed
 better in simulations.  In this method, missing data matches perfectly with missing data,
 while missing data matches partially with non-missing data.\cr
 
 \code{missingMethod=1} is retained for experimental purposes.  Here, missing data
 matches partially with missing and non-missing data. 
}
\value{
A distance/dissimilarity matrix of S3 class "amMatrix".
}
\references{
Please see the supplementary documentation for more information.  This is available as a vignette.  Click on the index link at the bottom of this page to find it.
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
myDataset <- amDataset(amExample1, missingCode="-99", indexColumn=1,
    metaDataColumn=2)

## Produce dissimilarity matrix
dissimMatrix <- amMatrix(myDataset)
}

}

