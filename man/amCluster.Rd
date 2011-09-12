\name{amCluster}
\alias{amCluster}
\alias{summary.amCluster}
\alias{amHTML.amCluster}
\alias{amCSV.amCluster}
\title{
Clustering of multilocus genotypes
}
\description{
Functions to perform clustering of multilocus genotypes in order to identify unique consensus and singleton genotypes, and review the output
of the analysis in formatted text, HTML, or CSV.  These functions are usually called by \code{\link{amUnique}}.  This interface remains to enable a better understanding of how \code{\link{amUnique}} operates.  For more information see example.

There are three steps to this analysis: (1) Identify the dissimilarity between pairs of genotypes using a metric
which takes missing data into account; (2) Cluster this dissimilarity matrix using a standard hierarchical agglomerative clustering approach; and (3) Use a dynamic tree cutting approach to identify clusters.   
}
\usage{
amCluster(amDatasetFocal, runUntilSingletons = TRUE, cutHeight = 0.3,
    missingMethod = 2, consensusMethod = 1, clusterMethod = "complete")

amHTML.amCluster(x, htmlFile = NULL, htmlCSS = amCSSForHTML())

amCSV.amCluster(x, csvFile)

\method{summary}{amCluster}(object, html = NULL, csv = NULL, ...)
}
\arguments{
  \item{amDatasetFocal}{
  An amDataset object containing genotypes to cluster
}
  \item{runUntilSingletons}{
  When \code{TRUE} the analysis runs recursively with the unique individuals determined in one analysis feeding into the next until no more clusters are formed.
  This is used when the goal is to thin a dataset to its unique individuals.  For more control over this process this can be done manually with \code{FALSE}.
  See details and examples.  
}
  \item{cutHeight}{
  Sets the tree cutting height using the hybrid method in the \code{\link{dynamicTreeCut}} package.  See details and \code{\link{cutreeHybrid}} for more information.
}
  \item{missingMethod}{
  The method used to determine the similarity of multilocus genotypes when data is
  missing.  The default (=2) is preferable in all cases.  Please see \code{\link{amMatrix}}.
}
  \item{consensusMethod}{
  The method (an integer) used to determine the consensus multilocus genotype from a cluster of multilocus genotypes.  See details.
}
  \item{clusterMethod}{
  The method used by \code{\link{hclust}} for clustering.  Only the default \code{clusterMethod = "complete"} performs acceptably in simulations.  This option remains for experimental reasons.
  
  
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
  The path to a CSV file to create containing only the unique genotypes determined in the clustering.

}
\item{...}{
Additional arguments to summary
}  
}
\details{
  
   Selecting an appropriate \code{cutHeight} parameter (also known as the d-hat criterion) is essential.  Typically this function is called from \code{amUnique}, and the conversion between alleleMismatch (m-hat) and cutHeight (d-hat) will be done automatically.  Selecting an appropriate value for alleleMismatch (m-hat) can be done using \code{amUniqueProfile}.  See the supplementary documentation for an explanation of how these parameters rae related.\cr

  \code{runUntilSingletons=TRUE} provides an efficient and reliable way to determine the unique individuals in a dataset if the dataset meets certain
  criteria.  To understand how the clustering is thinning the dataset run this recursion manually using
  \code{runUntilSingletons=FALSE}.  An example is given below.\cr

  \code{cutHeight} in practice gives the amount of dissimilarity (using the metric described in \code{\link{amMatrix}}) required for two multilocus
  genotypes to be declared different (also known as d-hat).  

  The default setting for \code{consensusMethod} performs well.
  \tabular{ll}{
  \code{consensusMethod} \tab  \cr
  \code{1} \tab Genotype with max similarity to others in the cluster is consensus (DEFAULT)\cr
  \code{2} \tab Genotype with max similarity to others in the cluster is consensus\cr
           \tab then interpolate missing alleles using mode non-missing allele in each column\cr
  \code{3} \tab Genotype with min missing data is consensus\cr
  \code{4} \tab Genotype with min missing data is consensus\cr
           \tab then interpolate missing alleles using mode non-missing allele in each column\cr
    }


  
}
\value{
An amCluster object.\cr
Or side effects:  analysis summary written to an HTML file or to the console, or unique genotypes written to a CSV file.
}

\note{
  There is an additional side effect of \code{html=TRUE} (or of \code{htmlFile=NULL}).  If required, there is a clean up of the
  operating system temporary directory where AlleleMatch temporary HTML files are stored.  Files that match the pattern am*.htm
  and are older 24 hours are deleted from this temporary directory.
}
\references{
Please see the supplementary documentation for more information.  This is available as a vignette.  Click on the index link at the bottom of this page to find it.
}
\author{
Paul Galpern (\email{pgalpern@gmail.com})
}
\seealso{
\code{\link{amDataset}}, \code{\link{amMatrix}}, \code{\link{amPairwise}}, \code{\link{amUnique}}, \code{\link{amUniqueProfile}}
}
\examples{
data("amExample5")

## Produce amDataset object
myDataset <- amDataset(amExample5, missingCode="-99", indexColumn=1,
    metaDataColumn=2, ignoreColumn="gender")

## Usage
myCluster <- amCluster(myDataset, cutHeight=0.2)

## Display analysis as HTML in default browser
summary(myCluster, html=TRUE)

## Save analysis to HTML file
summary(myCluster, html="myCluster.htm")

## Display analysis as formatted text on the console
summary(myCluster)

## Save unique genotypes only to a CSV file
summary(myCluster, csv="myCluster.csv")

## Demonstration of how amCluster operates
## Manual control over the recursion in amCluster()
summary(myCluster1 <- amCluster(myDataset, runUntilSingletons=FALSE,
    cutHeight=0.2), html=TRUE)
summary(myCluster2 <- amCluster(myCluster1$unique, runUntilSingletons=FALSE,
    cutHeight=0.2), html=TRUE)
summary(myCluster3 <- amCluster(myCluster2$unique, runUntilSingletons=FALSE,
    cutHeight=0.2), html=TRUE)
summary(myCluster4 <- amCluster(myCluster3$unique, runUntilSingletons=FALSE,
    cutHeight=0.2), html=TRUE)
## No more clusters, therefore stop.




}
