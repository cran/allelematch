\name{amPairwise}
\alias{amPairwise}
\alias{summary.amPairwise}
\alias{amHTML.amPairwise}
\alias{amCSV.amPairwise}
\title{
Pairwise matching of multilocus genotypes
}
\description{
Functions to perform a pairwise matching analysis of a multilocus genotype dataset, and review the output
in formatted text or HTML.  For each genotype in the focal dataset all genotypes in the comparison genotype
are returned that match at or above a threshold matching score.  The matching score is also known as the s-hat
criterion (see the supplementary documentation).  This is determined using \code{\link{amMatrix}}.
}
\usage{
amPairwise(amDatasetFocal, amDatasetComparison = amDatasetFocal,
    alleleMismatch=NULL, matchThreshold = NULL, missingMethod = 2)

amHTML.amPairwise(x, htmlFile = NULL, htmlCSS = amCSSForHTML())

amCSV.amPairwise(x, csvFile)

\method{summary}{amPairwise}(object, html = NULL, csv = NULL, ...)
}
\arguments{
  \item{amDatasetFocal}{
  An amDataset object containing focal genotypes.
}
  \item{amDatasetComparison}{
  Optional.  An amDataset object containing comparison genotypes.  If not supplied, the focal dataset
  is also the comparison dataset (i.e. all members of the focal dataset are compared against themselves).
}
\item{alleleMismatch}{
  Maximum number of mismatching alleles which will be tolerated when identifying individuals.  Also known as m-hat parameter.  If given, \code{matchThreshold} should be omitted.  
}
  \item{matchThreshold}{
  Return comparison genotypes that match with the focal genotype at or above this score or similarity.  Also known as s-hat parameter.
}
  \item{missingMethod}{
  The method used to determine the similarity of multilocus genotypes when data is
  missing.  The default (=2) is preferable in all cases.  Please see \code{\link{amMatrix}}.
  
  
}
 \item{object, x}{
 An amPairwise object.
}
  \item{htmlFile}{
  The path to an HTML file to create.  If \code{htmlFile=NULL} a file is created in the operating system
  temporary directory and is then opened in the default browser.
}
  \item{htmlCSS}{
  A string containing a valid cascading style sheet.  A default style sheet is provided in
  \code{\link{amCSSForHTML}}. See \code{\link{amCSSForHTML}} for details of how to tweak
  this CSS.
}
  \item{html}{
  If \code{html=NULL} or \code{html=FALSE} formatted textual output is displayed on the console.\cr
  If \code{html=TRUE} the summary method produces and loads an HTML file in the default browser.\cr
  \code{html}  can also contain a path to a file where HTML output will be written.
}
\item{csvFile, csv}{
  The path to a CSV file to create giving a data frame representation of the pairwise matching results.
}
\item{...}{
Additional arguments to summary
}  
}
\details{
Pairwise matching of genotypes is a useful means to assess data quality and inspect for genotyping errors.\cr

The \code{matchThreshold} represents the similarity between two multilocus genotypes, and can be thought of
as a percentage similarity (or a Hamming's distance between two vectors) that has been
corrected where missing data is present, such that missing data represents neither a match nor a mismatch, but
a "partial" match.  Please see \code{\link{amMatrix}} for more discussion of this metric.
}
\note{
  As \code{matchThreshold} is lowered the size of the output increases rapidly.  Typically analyses will not
  be very useful or manageable with thresholds below 0.7.\cr
  
  There is an additional side effect of \code{html=TRUE} (or of \code{htmlFile=NULL}).  If required, there is a clean up of the
  operating system temporary directory where AlleleMatch temporary HTML files are stored.  Files that match the pattern am*.htm
  and are older 24 hours are deleted from this temporary directory.
}
\value{
An amPairwise object.\cr
Or side effects:  analysis summary written to an HTML file or to the console.
}
\references{
Please see the supplementary documentation for more information.  This is available as a vignette.  Click on the index link at the bottom of this page to find it.
}
\author{
Paul Galpern (\email{pgalpern@gmail.com})
}
\seealso{
\code{\link{amDataset}}, \code{\link{amMatrix}}, \code{\link{amUnique}}
}
\examples{

\dontrun{

data("amExample5")

## Produce amDataset object
myDataset <- amDataset(amExample5, missingCode="-99", indexColumn=1,
    metaDataColumn=2, ignoreColumn="gender")

## Typical usage
myPairwise <- amPairwise(myDataset, alleleMismatch=2)

## Display analysis as HTML in default browser
summary(myPairwise, html=TRUE)

## Save analysis to HTML file
summary(myPairwise, html="myPairwise.htm")

## Save analysis to CSV file
summary(myPairwise, csv="myPairwise.csv")

## Display analysis as formatted text on the console
summary(myPairwise)

## Compare one dataset against a second
## Both must have same number of allele columns
## Here we create two datasets artificially from one for illustration purposes
myDatasetA <- amDataset(amExample5[sample(nrow(amExample5))[1:25], ],
    missingCode="-99", indexColumn=1, ignoreColumn=2)
myDatasetB <- amDataset(amExample5[sample(nrow(amExample5))[1:100], ],
    missingCode="-99", indexColumn=1, ignoreColumn=2)
myPairwise2 <- amPairwise(myDatasetA, myDatasetB, alleleMismatch=3)
summary(myPairwise2, html=TRUE)

}

}
