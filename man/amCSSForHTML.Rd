\name{amCSSForHTML}
\alias{amCSSForHTML}
\title{
Produce cascading style sheet (CSS) for HTML
}
\description{
This returns a string containing the CSS code for embedding in HTML output by \code{\link{amHTML.amPairwise}} and
\code{\link{amHTML.amUnique}}.  
}
\usage{
amCSSForHTML()
}

\details{
This function is used internally.  It can also be used as a basis to tweak the CSS
code if different output formatting and colour-coding are desired.  Please see examples.
}
\value{
A string containing a complete cascading style sheet
}
\references{
Please see the supplementary documentation for more information.  This is available as a vignette.  Click on the index link at the bottom of this page to find it.
}
\author{
Paul Galpern (\email{pgalpern@gmail.com})
}
\seealso{
\code{\link{amHTML.amPairwise}}, \code{\link{amHTML.amUnique}}
}
\examples{

data("amExample5")

## Produce CSS file for editing
cat(amCSSForHTML(), file="myCSS.css")

## Edit myCSS.css manually using text editor, then reimport as follows
myCSS <- paste(readLines("myCSS.css"), collapse="\n")

## Perform an allelematch unique analysis
myUnique <- amUnique(amDataset(amExample5, missingCode="-99", indexColumn=1,
    metaDataColumn=2, ignoreColumn="gender"),
alleleMismatch=3)

## Produce HTML output using tweaked CSS
amHTML.amUnique(myUnique, html="myUnique.htm", htmlCSS=myCSS)

}